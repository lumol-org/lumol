import subprocess
import time

COMMITS = [9, 8]


def run_benchmarks(out_name):
    for _ in range(3):
        with open('bench_results/%s' % out_name, 'w') as out:
            subprocess.call('cargo bench', shell=True, stdout=out)


def main():
    for i, n in enumerate(COMMITS):
        run_benchmarks('commit_u_%s' % n)
        if i > 0:
            for _ in range(COMMITS[i - 1] - n):
                subprocess.call('git checkout HEAD~', shell=True)

    subprocess.call('git checkout boost-ewald-forces', shell=True)


if __name__ == '__main__':
    main()
