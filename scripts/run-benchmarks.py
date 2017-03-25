import subprocess


def check_if_benches_should_run():
    EXPR = '__TRAVIS_RUN_BENCHES='
    out = subprocess.check_output('git log -n 1', shell=True).decode('utf-8')

    for line in out.split('\n'):
        # If the EXPR appears in the line, return the given value
        try:
            idx = line.index(EXPR)
            return int(line[idx + len(EXPR):])
        except ValueError:
            continue

    return None


def run_warmup():
    print('=================== Warming up ==============================')
    for _ in range(3):
        subprocess.check_output('cargo bench', shell=True)


def run_benches():
    subprocess.call('cargo bench', shell=True)


def run_benches_for_n_commits(n_commits):
    # Store the name of the current branch
    current_branch_name = subprocess \
        .check_output('git rev-parse --abbrev-ref HEAD', shell=True) \
        .decode('utf-8').rstrip()

    for k in range(n_commits):
        print('=================== Benching commit %s ==============================' % k)
        run_benches()
        subprocess.call('git checkout HEAD~', shell=True)
        print('=================== Done ==============================')

    # Go back to the commit handled by the CI
    subprocess.call('git checkout %s' % current_branch_name, shell=True)


def main():
    n_commits = check_if_benches_should_run()
    if n_commits is not None:
        run_warmup()
        run_benches_for_n_commits(n_commits)
    else:
        print("No benchmarks were run.")


if __name__ == '__main__':
    main()
