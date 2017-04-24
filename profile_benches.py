import glob
import os
import subprocess
import io

# files = glob.glob('benches/*.rs')
# files = [n.split('/')[-1].split('.')[0] for n in files]

files = ['water', 'nacl']
files.sort()


def get_executable_files(name):
    query = 'target/release/deps/%s-*' % name
    return glob.glob(query)


cmd = ['cargo', 'bench', '--no-run']
subprocess.call(cmd)
# if False:
#     for n in files:
#         print(n)
#         # if len(result) == 0:
#         #     continue
#         for name in get_executable_files(n):
#             os.remove(name)
#         cmd = ['cargo', 'build', '--bench', n, '--release']
#         # print(' '.join(cmd))
#         subprocess.call(cmd)
# os.mkdir('flamegraphs')
for n in files:
    # cmd = ['cargo', 'bench', '--bench', n, 'forces']
    # subprocess.call(cmd)
    exe_name = get_executable_files(n)[0]
    #
    # # exe_file = result[0]
    cmd = ['perf', 'record', '-g', './' + exe_name, 'virial_ewald']
    # for _ in range(3):
    with open('perf.data', 'w') as out:
        print(' '.join(cmd))
        subprocess.call(cmd, stdout=out)

    print(' '.join(['perf', 'script']))
    # out = io.StringIO()
    # p1 = subprocess.call(['perf', 'script', '--input', 'perf.data'], shell=True)
    subprocess.call('perf script --input=perf.data | stackcollapse-perf | flamegraph > flamegraphs/%s.html' % n,
                    shell=True)
    # for line in p1.stdout:
    #     print(line)
    # print(out)
    # print(' '.join(['perf', 'script']))
    # subprocess.call(['stackcollapse-perf'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    # with open('%s.html' % n, 'w') as out:
    #     subprocess.call(['flamegraph'], stdout=out)

    # break
    # subprocess.call(['perf', 'report'])
