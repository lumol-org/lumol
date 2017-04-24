import subprocess

for i in range(1, 4):
    subprocess.call('cargo benchcmp bench_results/commit_%s bench_results/commit_%s' % (i, i + 1),
                    shell=True)
