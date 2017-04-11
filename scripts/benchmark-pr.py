#!/usr/bin/env python
# coding=utf-8

import subprocess
import os
from contextlib import contextmanager

import click
import cpuinfo
import requests

BASE_REPO_URL = 'https://github.com/lumol-org/lumol.git'

INDIVIDUAL_BENCH_TEMPLATE = """
<details><summary>{sha} {title}</summary>
  <p>\n\n
\n\n
```bash\n
{body}\n\n
```\n\n

</p></details>
"""

COMPARISON_TEMPLATE = """
<details><summary>{sha} {title}</summary>
  <p>\n\n
\n\n
```bash\n
{body}\n\n
```\n\n

</p></details>
"""


def _get_environment_variable(name):
    var = os.environ.get(name, None)
    if var is None:
        msg = 'The environment variable {} must be set to access Github API'.format(name)
        raise Exception(msg)

    return var


def check_for_cargo_benchcmp():
    try:
        subprocess.check_output('cargo benchcmp --help', shell=True)
    except subprocess.CalledProcessError as e:
        raise Exception("cargo-benchcmp is not installed, please install with "
                        "`cargo install cargo-benchcmp`")


def check_for_environment_variables():
    """
    Used to check early that the required env vars are set.
    """
    _get_environment_variable('LUMOL_GH_USERNAME')
    _get_environment_variable('LUMOL_GH_TOKEN')


def request_api(endpoint, data=None):
    url = 'https://api.github.com/repos/lumol-org/lumol' + endpoint
    username = _get_environment_variable('LUMOL_GH_USERNAME')
    token = _get_environment_variable('LUMOL_GH_TOKEN')

    if data is None:
        r = requests.get(url, auth=(username, token))
    else:
        r = requests.post(url, auth=(username, token), json=data)

    return r.json()


def get_master_commit():
    cmd = 'git log --oneline _bot_remote/master -n 1'
    out = subprocess.check_output(cmd, shell=True).decode('utf-8')
    sha, _, title = out.rstrip().partition(' ')
    return sha, title


def get_commit_descriptions(n_commits=None):
    """
    Get hash and title of the `n_commits` latest commits on the PR.

    Also adds the commit at the HEAD of master in the end. If there
    are more than `n_commits` on the PR, only the ones on the PR
    are kept.

    :param n_commits:
    :return:
    """
    cmd = 'git log --oneline  _bot_remote/master.._bot_pr'.format(n_commits)
    out = subprocess.check_output(cmd, shell=True).decode('utf-8')

    descriptions = []
    for line in out.split('\n'):
        if line.rstrip() == '':
            continue
        sha, _, title = line.partition(' ')
        descriptions.append((sha, title))

    if n_commits is not None:
        descriptions = descriptions[:min(len(descriptions), n_commits)]

    descriptions.append(get_master_commit())
    return descriptions


def clean_repo():
    subprocess.call('git checkout master', shell=True)
    subprocess.call('git branch -D _bot_pr', shell=True)
    subprocess.call('git remote remove _bot_remote', shell=True)


@contextmanager
def setup_repo(pr_id):
    clean_repo()
    subprocess.call('git remote add _bot_remote {}'.format(BASE_REPO_URL), shell=True)
    subprocess.call('git fetch _bot_remote master'.format(pr_id), shell=True)
    subprocess.call('git fetch _bot_remote pull/{}/head:_bot_pr'.format(pr_id), shell=True)
    subprocess.call('git checkout _bot_pr', shell=True)
    try:
        yield
    finally:
        clean_repo()


class Benchmarker:
    def __init__(self, n_commits, output_dir):
        self.commit_descriptions = get_commit_descriptions(n_commits)
        self.output_dir = os.path.abspath(os.path.expanduser(output_dir))

    def run_warmup(self):
        print('=================== Warming up ==============================')
        for _ in range(3):
            subprocess.check_output('cargo bench', shell=True)

    def run_bench(self, sha):
        cmd = 'cargo bench > {}/{}.txt'.format(self.output_dir, sha)
        subprocess.call(cmd, shell=True)

    def run_all_benches(self):
        for sha, title in self.commit_descriptions:
            print('=================== Benchmarking commit {} =============================='.format(sha))
            subprocess.call('git checkout {}'.format(sha), shell=True)
            self.run_bench(sha)
            print('=================== Done ==============================')

    def compare_benches(self):
        comparisons = {}
        master_sha, _ = self.commit_descriptions[-1]
        for sha, title in self.commit_descriptions[:-1]:
            cmd = 'cargo benchcmp {dir}/{sha_1}.txt {dir}/{sha_2}.txt --threshold 2 --variance' \
                .format(dir=self.output_dir, sha_1=master_sha, sha_2=sha)
            out = subprocess.check_output(cmd, shell=True).decode('utf-8')
            comparisons[sha] = out

        return comparisons

    def comment_pr(self, pr_id):
        # Comparison benchmarks
        master_sha, master_title = self.commit_descriptions[-1]
        comment = '## Comparing to master ({})\nusing `--threshold 2, latest commit first`'.format(master_sha)

        comparisons = self.compare_benches()
        for sha, title in self.commit_descriptions[:-1]:
            compare = comparisons[sha]
            comment += COMPARISON_TEMPLATE.format(sha=sha, title=title, body=compare)

        # Individual benchmarks
        comment += '\n## Individual benchmarks\n'

        for k, (sha, title) in enumerate(self.commit_descriptions):
            with open('{}/{}.txt'.format(self.output_dir, sha)) as f:
                bench = f.read()
            comment += INDIVIDUAL_BENCH_TEMPLATE.format(sha=sha, title=title, body=bench)

        info = cpuinfo.get_cpu_info()
        if info is not None:
            comment += '\n<br>**CPU**: {}'.format(info['brand'])

        # Emit the request
        data = {
            'body': comment
        }
        request_api('/issues/{}/comments'.format(pr_id), data)
        print('Comment posted on PR.')


@click.command()
@click.argument('pr_id', type=click.INT)
@click.argument('output_dir', default='./target/benchmarks/')
@click.option('--n-commits', '-n', type=click.INT, help='maximum number of commits to benchmark')
def main(pr_id, output_dir, n_commits):
    """
    Run the benchmarks for multiple commits on a PR and compare to master.

    The benchmark results are saved in OUTPUT_DIR, and a comment
    with a summary will be automatically added to the PR.
    This script requires the environment variables LUMOL_GH_USERNAME
    and LUMOL_GH_TOKEN to contain respectively the Github username
    and a personal access token.
    """
    check_for_cargo_benchcmp()
    check_for_environment_variables()

    os.makedirs(output_dir, exist_ok=True)

    with setup_repo(pr_id):
        benchmarker = Benchmarker(n_commits, output_dir)
        benchmarker.run_warmup()
        benchmarker.run_all_benches()
        benchmarker.comment_pr(pr_id)


if __name__ == '__main__':
    main()
