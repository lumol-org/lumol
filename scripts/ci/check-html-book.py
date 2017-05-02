#!/usr/bin/env python
# coding=utf-8
import os
import sys
from bs4 import BeautifulSoup


ERRORS = 0


def error(message):
    global ERRORS
    ERRORS += 1
    print(message)


def get_all_targets(dir):
    targets = set()
    for (root, _, paths) in os.walk(dir):
            for path in paths:
                # Link to the page
                targets.add(os.path.relpath(os.path.join(root, path), dir))

                extension = path.split(".")[-1]
                if extension == "html":
                    # Link to a heading in the page
                    html = open(os.path.join(root, path)).read()
                    soup = BeautifulSoup(html, 'html.parser')
                    for heading in soup.find_all("a", class_="header"):
                        targets.add(heading["href"])
    return targets


def check_links(path, targets):
    def is_local_link(link):
        return (
            link and
            not link.startswith("http://") and
            not link.startswith("https://")
        )

    html = open(path).read()
    soup = BeautifulSoup(html, 'html.parser')
    for link in soup.find_all('a'):
        target = link["href"]
        if is_local_link(target):
            if target not in targets:
                error("broken link to '{}' in {}".format(target, path))


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("usage: {} path/to/html/book".format(sys.argv[0]))
        sys.exit(1)

    html_book = sys.argv[1]
    targets = get_all_targets(html_book)
    for (root, _, paths) in os.walk(html_book):
            for path in paths:
                if path == "print.html" and root == html_book:
                    # Ignore the print.html file, as it contains the same
                    # HTML content
                    continue
                extension = path.split(".")[-1]
                if extension == "html":
                    check_links(os.path.join(root, path), targets)

    if ERRORS != 0:
        print("------------------\n{} link errors in doc".format(ERRORS))
        sys.exit(1)
