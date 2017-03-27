import subprocess
import json
import os

INDIVIDUAL_BENCH_TEMPLATE = """
<details><summary>%s %s</summary>
  <p>\n\n
\n\n
```bash\n
%s\n\n
```\n\n

</p></details>
"""

COMPARISON_TEMPLATE = """
<details><summary>%s %s</summary>
  <p>\n\n
\n\n
```bash\n
%s\n\n
```\n\n

</p></details>
"""


def get_pull_request_id():
    pr_id = os.environ.get('TRAVIS_PULL_REQUEST', 'false')
    if pr_id == 'false':
        return None
    return int(pr_id)


def request_api(endpoint, data=None):
    url = 'https://api.github.com/repos/lumol-org/lumol' + endpoint
    username = os.environ['GITHUB_USER_NAME']
    token = os.environ['GITHUB_USER_TOKEN']

    cmd = "curl -sS %s  -u '%s:%s'" % (url, username, token)
    if data is not None:
        cmd += " --data '%s'" % json.dumps(data)

    out = subprocess.check_output(cmd, shell=True).decode('utf-8')
    return json.loads(out)


def check_if_benches_should_run():
    EXPR = '@bot benchmark'
    pr_id = get_pull_request_id()
    if pr_id is None:
        print("This is not a PR, the benchmarks are not run.")
        return None

    # Get the text of the PR in `body`
    out = request_api('/pulls/%s' % pr_id)

    clone_url = out['head']['repo']['clone_url']
    body = out['body']

    if EXPR not in body:
        print('No benchmarks were requested.')
        return None

    begin = body.index(EXPR)
    end = len(body)
    try:
        s = begin + len(EXPR) + 1
        end = body[s:].index(' ') + s
    except ValueError:
        pass

    n_commits = body[begin + len(EXPR):end]

    return int(n_commits), clone_url


def get_commit_descriptions(n_commits):
    cmd = 'git log --oneline | head -n %s' % n_commits
    out = subprocess.check_output(cmd, shell=True).decode('utf-8')
    descriptions = []
    for line in out.split('\n'):
        if line.rstrip() == '':
            continue
        h, _, title = line.partition(' ')
        descriptions.append((h, title))

    return descriptions


def run_warmup():
    print('=================== Warming up ==============================')
    for _ in range(3):
        subprocess.check_output('cargo bench', shell=True)


def run_bench(k):
    # The output files have to be outside of the directory for
    # us to be able to git checkout.
    subprocess.call('cargo bench > ../bench_%s.txt' % k, shell=True)


def run_benches_for_n_commits(n_commits):
    for k in range(n_commits):
        print('=================== Benching commit %s ==============================' % k)
        run_bench(k)
        subprocess.call('git checkout HEAD~', shell=True)
        print('=================== Done ==============================')


def cd_to_cloned_repo(url):
    subprocess.call('git clone %s cloned_repo' % url, shell=True)
    os.chdir('cloned_repo')
    subprocess.call('pwd', shell=True)
    subprocess.call('git checkout %s' % os.environ['TRAVIS_PULL_REQUEST_BRANCH'], shell=True)
    subprocess.call('git status', shell=True)


def comment_pull_request(n_commits):
    pr_id = get_pull_request_id()
    descriptions = get_commit_descriptions(n_commits)

    # Comparison benchmarks
    comment = '## Comparison benchmarks\n (comparisons with %s)' % descriptions[-1][0]

    for k, (h, title) in enumerate(descriptions[:-1]):
        cmd = 'cargo benchcmp ../bench_%s.txt ../bench_%s' % (n_commits - 1, k)
        out = subprocess.check_output(cmd, shell=True).decode('utf-8')
        comment += COMPARISON_TEMPLATE % (h, title, out)

    # Individual benchmarks
    comment += '\n## Individual benchmarks\n'

    for k, (h, title) in enumerate(descriptions):
        with open('../bench_%s.txt' % k) as f:
            bench = f.read()
        comment += INDIVIDUAL_BENCH_TEMPLATE % (h, title, bench)

    # Emit the request
    data = {
        'body': comment
    }
    request_api('/issues/%s/comments' % pr_id, data)


def main():
    result = check_if_benches_should_run()
    if result is None:
        return
    n_commits, clone_url = result
    cd_to_cloned_repo(clone_url)
    run_warmup()
    run_benches_for_n_commits(n_commits)
    comment_pull_request(n_commits)
    subprocess.call('rm ../bench_*.txt', shell=True)


if __name__ == '__main__':
    main()
