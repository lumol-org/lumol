#!/usr/bin/env python3
# -*- encoding: utf8 -*-
import json
import subprocess


def list_executables(*command):
    command = list(command)
    command.append("--no-run")
    command.append("--message-format=json")

    result = subprocess.run(command, check=True, stdout=subprocess.PIPE)
    # convert successive JSON documents to a list of JSON objects
    stdout = "[" + result.stdout.decode("utf8").replace("}\n{", "},{") + "]"

    executables = []
    for event in json.loads(stdout):
        if "profile" in event and event["profile"]["test"]:
            for path in event["filenames"]:
                if path.endswith(".dSYM"):
                    continue
                executables.append(path)

    return executables


def main():
    executables = []
    executables.extend(list_executables("cargo", "test", "--all", "--lib"))
    executables.extend(
        list_executables("cargo", "test", "--package=lumol-input", "--tests")
    )

    for executable in set(executables):
        if "lumol_tutorial" in executable:
            continue
        print(executable)


if __name__ == "__main__":
    main()
