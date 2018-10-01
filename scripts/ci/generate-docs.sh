#!/bin/bash -ex
# Build the docs and the user manual

# Build the doc
RUSTDOCFLAGS="--html-in-header doc/KaTeX.html" cargo doc --no-deps
make -C doc html

# Move it to the right place
rm -rf target/gh-pages/latest
mkdir -p target/gh-pages/latest
cd target/gh-pages/
cp -r ../../doc/build/html latest/book
cp ../../doc/index.html latest/index.html
mv ../doc/* latest/

cat <<EOF > index.html
<meta http-equiv=refresh content=0;url=latest/index.html>
EOF

touch .nojekyll
