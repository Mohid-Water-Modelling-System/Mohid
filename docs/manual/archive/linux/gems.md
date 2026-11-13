# Gems

**Rubygems** is the extension manager for [ruby](ruby.md). Much like cpan is the extension manager for [Perl](perl.md). **Gems** is more than a packaging manager: it's a library versioning manager too! [That is also it's achille's heel](http://blade.nagaokaut.ac.jp/cgi-bin/scat.rb/ruby/ruby-talk/132899).

## Install

Requires [ruby](ruby.md) to be installed prior.

Install [rubygems](http://rubygems.org/)

```
> wget http://rubyforge.org/frs/download.php/35283/rubygems-1.1.1.tgz
> tar -xvf rubygems-1.1.1.tgz
> cd rubygems-1.1.1
> ruby setup.rb
```

Install extra-packages for gems to work correctly

```
> sudo apt-get install ruby1.8-dev
```

To configure gems on your system you need to set the environment variable RUBYOPT. This will preload "require gems" to every ruby script you run.

```
>vim ~/.bashrc
.bashrc> export RUBYOPT=rubygems
```

## Popular gems

### xmpp4r

xmpp is the open protocol behind the [jabber](README.md#jabber#jabber) service.

```
> sudo apt-get install libopenssl-ruby
> sudo gem install xmpp4r -y
```

Testing

```
> irb
irb> require 'xmpp4r/client'
```

#### xmpp4r-simple

```
> sudo gem install xmpp4r-simple -y
```

### twitter

```
> sudo gem install twitter -y
```

## External references

- [Ruby](ruby.md),
- [gems homepage](http://rubygems.org/)
