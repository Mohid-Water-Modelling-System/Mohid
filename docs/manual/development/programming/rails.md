# Rails

Rails or [Ruby](../../archive/linux/ruby.md) on Rails is a web-development framework in [ruby](../../archive/linux/ruby.md) that uses the [model-view-controller (MVC)](http://en.wikipedia.org/wiki/Model-view-controller) pattern.

## Install

You need to have ruby and gems installed. Sqlite is also recommended.

```
> sudo yum install sqlite3 sqlite3-devel
> sudo gem install rails
```

## Your first rails app

```
> rails mywebapp
> cd mywebapp
> script/controller welcome hello
edit, then save, the mywebapp/app/views/welcome/hello.html.erb ...
> script/server
```

Open up in your browser *[http://localhost:3000/welcome/hello](http://localhost:3000/welcome/hello)*. That's it!

## See also

- [Ruby](../../archive/linux/ruby.md)
- [Python](../../archive/linux/python.md)
- [Perl](../../archive/linux/perl.md)

## External links

- [Rails tutorial](http://wiki.rubyonrails.org/rails/pages/Tutorial)
