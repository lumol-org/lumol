# Input files

Cymbalum uses the [YAML](http://yaml.org/) format for all the input files. This
format have significant whitespace, and is easy to read by humans.

If you are not familiar with the YAML syntax, [this
page](https://learnxinyminutes.com/docs/yaml/) will teach you the basics.

<!-- TODO: add an introduction to YAML -->

All the keys should be written in lower case, but string values are case
insensitive. You should always use CamelCase for string value, *i.e.* you should
prefer `LennardJones` to `lennardjones`.

```yaml
# This is OK and completely equivalent
bar: Foo
bar: FOO
bar: foo

# This is not
Bar: Foo
BAR: Foo
```
