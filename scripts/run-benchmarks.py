import subprocess
import json
import os


def check_if_benches_should_run():
    EXPR = '__TRAVIS_RUN_BENCHES='
    out = subprocess.check_output('git log -n 1', shell=True).decode('utf-8')
    print(os.environ.get('TRAVIS_COMMIT_MESSAGE', 'NO COMMIT MSG'))
    print(os.environ.get('TRAVIS_COMMIT', 'NO COMMIT'))
    print(os.environ.get('TRAVIS_PULL_REQUEST_BRANCH', 'NO COMMIT'))

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

    os.mkdir('bench_results')

    for k in range(n_commits):
        print('=================== Benching commit %s ==============================' % k)
        run_benches()
        subprocess.call('git checkout HEAD~', shell=True)
        print('=================== Done ==============================')

    # Go back to the commit handled by the CI
    subprocess.call('git checkout %s' % current_branch_name, shell=True)


def get_pull_request_id():
    pr_id = os.environ.get('TRAVIS_PULL_REQUEST', 'false')
    if pr_id == 'false':
        return None
    return int(pr_id)


def cd_to_cloned_repo():
    url = 'git@github.com:lumol-org/lumol.git'
    subprocess.call('git clone %s cloned_repo' % url)
    os.chdir('cloned_repo')
    subprocess.call('git status')


def comment_pull_request(comment):
    pr_id = get_pull_request_id()
    print(pr_id)
    username = os.environ.get('GITHUB_USER_NAME', None)
    token = os.environ.get('GITHUB_USER_TOKEN', None)

    if username is None:
        print('GITHUB_USER_NAME is not set, unable to comment PR.')

    if token is None:
        print('GITHUB_USER_TOKEN is not set, unable to comment PR.')

    data = {
        'body': comment
    }
    url = 'https://api.github.com/repos/lumol-org/lumol/issues/%s/comments' % pr_id
    cmd = "curl -i %s --data '%s' -u '%s:%s'" % (url, json.dumps(data), username, token)
    print("Calling %s" % cmd)
    subprocess.call(cmd, shell=True)


def main():
    n_commits = check_if_benches_should_run()
    if n_commits is not None:
        if get_pull_request_id() is None:
            print("This is not a PR, the benchmarks are not ran.")
        cd_to_cloned_repo()
        comment_pull_request("Test comment from Travis CI")
        run_warmup()
        run_benches_for_n_commits(n_commits)
        os.chdir('../')
    else:
        print("No benchmarks were run.")


if __name__ == '__main__':
    main()
