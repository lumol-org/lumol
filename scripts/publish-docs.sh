#!/bin/bash -ex
# Build the docs and the user manual, and push them to github pages.

# Build the doc
cd src/core
cargo doc
cd ../input
cargo doc
cd ../..
mdbook build doc

# Move it to the right place
cd target/website
git checkout gh-pages

rm -rf latest
mkdir latest
mv ../../doc/book latest/book
cp ../../doc/index.html latest/index.html
mv ../doc/* latest/

cat <<EOF > index.html
<meta http-equiv=refresh content=0;url=latest/>
EOF

git config user.name "Travis CI"
git config user.email "luthaf@luthaf.fr"
git config push.default simple
git add --all .

# Skip push if there is no change
if git diff --cached --exit-code --quiet; then
    echo "No changes to the output on this push; exiting."
    exit 0
fi

git commit -m "[AUTO-COMMIT] Documentation update"
git push --force
