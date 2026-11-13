# Ruby

Ruby is (said to be) a great language. People say that programming is fun again. Everything in Ruby is an object. Strings, floats, integers, regexps, arrays are all objects, and all have their methods and fields.

## Starters

1. [Ruby in 20 minute](http://www.ruby-lang.org/en/documentation/quickstart/)

Install ruby

```
> sudo apt-get install ruby irb rdoc
```

## Features

Here we list the interesting features that make ruby so interesting, besides its interesting programming features. These are the extra tools that make working with **ruby** appealing and easy.

1. interactive shell (irb),
2. extension packages manager (gems),
3. documentation generator from in-line comments (rdoc),
4. unit-testing framework (utest),
5. web-apps development framework (rails),

[Python](python.md) also possesses all of these features built-in (ipython, install, pydoc, unit-tests and django).

### Irb

The most interesting thing in shell scripting ([bash](bash.md), [dos](../../building/windows/batch.md)) is that you can play right from the shell command line. This is what makes interpreted languages so interesting in the first place. Thus, like [matlab](matlab.md) or [python](python.md) (however, I'm not aware of one for [perl](perl.md)), ruby has an interactive shell called irb. Interactive ruby allows to experiment and debug interactively in ruby. Similar (but not so powerful) to [ipython](python.md).

```
> irb
irb(main):001:0> puts "Hello world!"
Hello world!
=> nil
```

The most interesting method to call in irb is *public\_methods* which lists all the object's available methods:

```
irb(main):001:0> puts 11.public_methods.sort.join(" | ")
% | & | * | ** | + | +@ | - | -@ | / | < | << | <= | <=>
| == | === | =~ | > | >= | >> | [] | ^ |  __id__ | __send__
| abs | between? | ceil | chr | class | clone | coerce
| denominator | display |  div | divmod | downto | dup | eql?
| equal? | extend | floor | freeze | frozen? | gcd | gcdlcm
|  gem | hash | id | id2name | inspect | instance_eval | instance_of?
| instance_variable_defined? |  instance_variable_get | instance_variable_set
| instance_variables | integer? | is_a? | kind_of? |  lcm | method | methods
| modulo | next | nil? | nonzero? | numerator | object_id | power! | prec
|  prec_f | prec_i | pretty_inspect | pretty_print | pretty_print_cycle
| pretty_print_inspect |  pretty_print_instance_variables | private_methods
| protected_methods | public_methods | quo | rdiv | remainder | require
| require_gem | respond_to? | round | rpower | send | singleton_method_added
| singleton_methods | size | step | succ | taguri | taguri= | taint | tainted?
| times | to_a | to_bn | to_f | to_i | to_int | to_r | to_s | to_sym | to_yaml
| to_yaml_properties | to_yaml_style | truncate | type | untaint | upto | zero? | | | ~
=> nil
```

### [Gems](gems.md)

Requires [rubygems](http://rubygems.org/) to be installed. [Rubygems](gems.md) is the extension manager for ruby. Much like cpan is the extension manager for [Perl](perl.md). [Gems](gems.md) is more than a packaging manager: it's a library versioning manager too! [That is also it's achille's heel](http://blade.nagaokaut.ac.jp/cgi-bin/scat.rb/ruby/ruby-talk/132899).

### Rdoc

Rdoc is the automatic html documentation generator for ruby scripts. Much like javadoc is for java. It requires rdoc to be installed.

```
> rdoc --main mywork/mywork.rb
```

### Unit testing

Unit tests are XP's method of choice to build safe long-lasting software. Still looking to see how these work...

### Rails

Hopefully we'll get there soon enough...

### Ruby2exe

[rubyscript2exe](http://www.erikveen.dds.nl/rubyscript2exe/index.html#6.0.0) allows to transform ruby scripts into full-bloated executable programs.

```
> gem install rubyscript2exe
```

## Sample code

### loop

```
object = Array.new(["one", "two", "three", "four"])
object.each do |a|
  puts a + "\n"
end
```

### array

```
object = ["one", "two", "three", "four"]
target = []
object.each { |a|
  target.push(a)
}
puts target
```

### i/o

```
if File.exist? "sample.txt"
  File.open("sample.txt", "r") { |fs|
    while line = fs.gets
      puts "Line read was '#{line.chomp}'.\n"
    end
    fs.close
  }
end
```

```
File.open("sample.txt", "w") { |fw|
  fw.puts "Hi, this a new line.\n"
  fw.close
}
```

### regexp

```
pattern = Regexp.new("[a-z]+", "i")
#or
pattern = /[a-z]+/i
if File.exist? "sample.txt"
  File.open("sample.txt", "r") { |fs|
    while line = fs.gets
      if line.match(pattern)
        puts "Line read was '#{line.chomp}'.\n"
      end
    end
    fs.close
  }
end
```

### method

```
def whatisit?(arg1, arg2="optional")
  #Ducktyping: Does it walk like a duck? Does quacks like a duck?
  if arg1.responds_to? :capitalize
    puts "Then '#{arg1.capitalize}' must be a string.\n"
  else
    if arg1.responds_to? :modulo && arg1.responds_to? :numerator
      puts "Then '#{arg1.to_s}' must be an integer.\n"
    else
      if arg1.responds_to? :modulo && arg1.responds_to? :truncate
        puts "Then '#{arg1.to_s}' must be a float.\n"
      end
    end
  end
  puts "Input '#{arg2}' was completely optional.\n"
end
```

### class

```
#File checktype.rb
class Checktype

  #Constructor
  def initialize(arg)
    #class Checktype field
    @variable = arg
  end

  #class Checktype method
  def isWhat?(arg=@variable)
    #Ducktyping: Does it walk like a duck? Does it quack like a duck?
    if arg.responds_to? :capitalize
      puts "Then '#{arg.capitalize}' must be a string.\n"
    else
      if arg.responds_to? :modulo && arg.responds_to? :numerator
        puts "Then '#{arg.to_s}' must be an integer.\n"
      else
        if arg.responds_to? :modulo && arg.responds_to? :truncate
          puts "Then '#{arg.to_s}' must be a float.\n"
        end
      end
    end
  end

end
```

```
require 'checktype'
somevar = Checktype(23)
somevar.isWhat?
#=> '23' must be an integer
Checktype.isWhat? "Hi"
#=> 'Hi' must be a string
```

### Complex examples

- [RootTransform classes](../../development/programming/ruby-examples.md)

## External References

- [ruby-lang](http://www.ruby-lang.org/en/),
- [Programming Ruby](http://www.ruby-doc.org/docs/ProgrammingRuby/)
