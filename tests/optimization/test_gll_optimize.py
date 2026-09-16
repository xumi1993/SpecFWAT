#!/usr/bin/env python3
"""MPI smoke tests using a CMake Makefiles build and GNU-compatible linker.

Run with the build's compiler/MPI modules loaded:
    python3 tests/optimization/test_gll_optimize.py --build build

Only configuration/mesh inputs and expensive forward solves are replaced.
The production driver, parser, optimizer, misfit reader and binary I/O execute.
All generated models stay in a temporary directory; no case data are modified.
"""

import argparse
import array
import math
import os
from pathlib import Path
import shlex
import struct
import subprocess
import tempfile


def run(args, **kwargs):
    return subprocess.run(args, check=True, text=True, **kwargs)


def read_vector(path):
    data = path.read_bytes()
    size = struct.unpack('=i', data[:4])[0]
    assert size == len(data) - 8 == struct.unpack('=i', data[-4:])[0], path
    values = array.array('f')
    values.frombytes(data[4:-4])
    assert all(math.isfinite(v) for v in values), path
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', type=Path, default=Path('build'))
    parser.add_argument('--mpiexec', default='mpirun')
    args = parser.parse_args()
    build = args.build.resolve()
    source = Path(__file__).resolve().parent
    link_dir = build / 'src/optimization'
    run(['cmake', '--build', str(build), '--target', 'fwat_optimize', '-j', '4'])
    link = shlex.split((link_dir / 'CMakeFiles/fwat_optimize.dir/link.txt').read_text())
    with tempfile.TemporaryDirectory(prefix='fwat-optimize-') as temporary:
        work = Path(temporary)
        fixture = work / 'fixture.o'
        run([link[0], '-I' + str(build / 'modules'), '-c', str(source / 'gll_fixture.f90'),
             '-o', str(fixture)], cwd=work)
        symbols = subprocess.check_output(['nm', '--defined-only', str(fixture)], text=True)

        def symbol(module, procedure):
            if 'gll_fixture_mp_' in symbols:
                return module + '_mp_' + procedure + '_'
            return '__' + module + '_MOD_' + procedure

        replacements = {
            'start_mpi': 'init_mpi_',
            'read_config': symbol('input_params', 'read_fwat_parameter_file'),
            'read_parameters': 'read_parameter_file_',
            'read_mesh': symbol('kernel_io', 'read_mesh_databases_for_init'),
            'check_resolution': 'check_mesh_resolution_',
            'read_events': symbol('input_params', 'acqui_read_source_set'),
            'free_events': symbol('input_params', 'acqui_finalize'),
            'generate_databases': symbol('generate_databases_subs', 'generate_databases_fwat'),
            'forward': symbol('line_search', 'forward_for_simu_type'),
        }
        rename = ['objcopy']
        for name, target in replacements.items():
            rename += ['--redefine-sym', symbol('gll_fixture', name) + '=__wrap_' + target]
        run(rename + [str(fixture)])
        executable = work / 'optimize-test'
        link[link.index('-o') + 1] = str(executable)
        link += [str(fixture)] + ['-Wl,--wrap=' + name for name in replacements.values()]
        run(link, cwd=link_dir)

        count = 0

        def check(case, options, message=None):
            nonlocal count
            count += 1
            directory = work / str(count)
            directory.mkdir()
            (directory / 'OUTPUT_FILES').mkdir()
            env = dict(os.environ, GLL_TEST_CASE=case, OMP_NUM_THREADS='1')
            result = subprocess.run([args.mpiexec, '-np', '2', str(executable), *options],
                                    cwd=directory, env=env, text=True, capture_output=True, timeout=60)
            output = result.stdout + result.stderr
            if message is not None:
                assert result.returncode != 0 and message in output, output
                assert not (directory / 'optimize').exists(), 'Invalid input wrote models'
            else:
                assert result.returncode == 0, output
            print('PASS:', case, ' '.join(options))
            return directory, output

        for case in ('sd', 'external', 'precond', 'lbfgs', 'zero', 'linesearch'):
            current = 'M01' if case == 'lbfgs' else 'M00'
            following = 'M02' if case == 'lbfgs' else 'M01'
            directory, _ = check(case, ['-m', current, '-g'])
            for rank in range(2):
                old, new, direction = [], [], []
                for name in ('vp', 'vs', 'rho'):
                    filename = f'proc{rank:06d}_{name}.bin'
                    old.append(read_vector(directory / 'optimize' / ('model_' + current) / filename))
                    new.append(read_vector(directory / 'DATABASES_MPI' / filename))
                    if case != 'zero':
                        archived = directory / 'optimize' / ('model_' + following) / filename
                        assert read_vector(archived) == new[-1]
                        direction.append(read_vector(directory / 'optimize' / ('DIRECTION_' + current) / filename))
                if case == 'zero':
                    assert old == new
                    assert not (directory / 'optimize' / ('model_' + following)).exists()
                    continue
                step = 0.05 if case == 'linesearch' else 0.1
                for point in range(len(old[0])):
                    expected = [old[p][point] * math.exp(step * direction[p][point]) for p in range(3)]
                    expected[0] = max(1.3 * expected[1], min(2.5 * expected[1], expected[0]))
                    for p in range(3):
                        assert math.isclose(new[p][point], expected[p], rel_tol=2e-6)
                if case not in ('lbfgs',):
                    factor = (2.0 if rank == 0 else 0.5) if case == 'precond' else 1.0
                    for p in range(3):
                        expected_direction = -(p + 1) * (rank + 1) * factor / 6.0
                        assert math.isclose(direction[p][0], expected_direction, rel_tol=2e-6)

        check('sd', ['--use-gll', '--model', 'M00'])
        _, output = check('sd', ['-m', 'M00'])
        assert 'regular-grid mode selected' in output
        _, output = check('sd', ['--help'])
        assert '--use-gll' in output
        check('sd', ['-g'], 'Model name not set')
        check('sd', ['-m', '-g'], 'Model name not set')
        check('sd', ['-m', 'M00', '--unknown'], 'Unknown option')
        check('sd', ['-m', 'M00', '-m', 'M01'], 'Repeated model option')
        check('sd', ['-m', 'M99', '-g'], 'between ITER_START and M98')
        check('sd', ['-m', 'oops', '-g'], 'form M00 through M98')
        check('joint', ['-m', 'M00', '-g'], 'exactly one POSTPROC.INV_TYPE')
        check('bad_method', ['-m', 'M00', '-g'], 'OPT_METHOD: 1 (SD) or 2')
        check('bad_mpi', ['-m', 'M00', '-g'], 'exactly NPROC ranks')
        check('external', ['-m', 'M01', '-g'], 'Set MODEL = gll')
        print(f'All {count} MPI smoke cases passed.')


if __name__ == '__main__':
    main()
