# Vim

Here are some useful references for using **Vim**, the powerful customizable wrist-friendly text-editor:

## Basic editing

**

$$
ESC
$$

** to return in normal mode.

**i** key to insert text.

**h, j, k, l** keys for moving around the text.

**x** key to delete a character.

**rx** keys replaces current character with the x character.

**u** key to undo an editing operation.

**U** key to undo all editing operations performed to the current line.

**

$$
Ctrl
$$

-r** keys to redo (undo the undo).

## Motion keys

**$** key to move to the end of line.

**0** key to move to the beginning of line.

**w** key to move to the first character of the next word.

**e** key to move to the last character of the current(or next) word.

## Operators

Operators are useful to perform operations a repeated number of times. Syntax is

```
Op [count] motion
```

wher Op is the operator, 

$$
count
$$

 is the number of times, and motion is the scope of the operation.

### Delete

**d** key.

**dd** deletes the whole line.

### Put

Whenever a delete operation is performed, the deleted characters are stored in the put variable.

**p** key.

```
example: 'dd' followed by 'p' performs a cut/paste one-line operation
```

### Change

**c** key.

## Cursor location

**

$$
CTRL-g
$$

** shows the cursor current position in text.

**gg** Goes to the beginning of text.

**G** Goes to the end of text.

**

$$
line
$$

G** goes to the given line.

## Search

**/text** keys finds the first occurence of *text* in the buffer from the cursor.

**?text** keys finds the first backward occurence of *text* in the buffer from the cursor.

**n** key finds the next occurence of *text*.

**N** key finds the previous occurence of *text*.

**

$$
Ctrl-o
$$

**

**

$$
Ctrl-i
$$

**

### Replace

**s/old text/new text** replaces first occurence in line.

**s/old text/new text/g** replaces all in line.

**#,#s/old text/new text/g** replaces all between lines # and #.

**%s/old text/new text/g** replaces all in buffer.

**%s/old text/new text/gc** replaces all in buffer but prompts before each replacement.

### Replace in several files

You need two commands. The first uses wildcards to select the group of files to be added in the "args". The second applies a search-and-replace-in-all-text to the "args".

```
:args ./*Base1/*.F90
:argdo %s@FOO@BAR@ge | update
```

## Add ftplugin

Copy the plugin into the `$VIM/ftplugin` directory then type at the vim prompt

```
:filetype plugin on
```

## Configure VIM

Edit the `/etc/vimrc` file or the `~/.vimrc` in linux or the `$VIM/_vimrc` file in windows.

### Space indentation

This will correctly set <tabs> to 4 spaces.

```
set shiftwidth=4        " defines the indentation characters
set tabstop=4           " defines the indentation characters
set expandtab           " defines the indentation characters
```

### Syntax highlighting

```
:syntax on
```

#### Adding a new language in syntax highlighting

```
> cp [lang].vim ~/.vim/syntax/[lang].vim
> locate filetype.vim
> sudo vim $PATHTOFILETYPE/filetype.vim
filetype.vim> au BufNewFile,BufRead *.[lang_suf]       setf [lang]
```

### Code completion (intellisense) for C++

Installing auto-completion for C++ in Vim

```
> sudo apt-get install vim exuberant-ctags
> wget [omnicppcomplete]
> tar -xvf [omnicppcomplete] ~/.vim
> vim ~/.vimrc
vim> :helptags ~/.vim/doc
vim> [append in ~/.vimrc]
" --- OmniCppComplete ---
" -- required --
set nocp " non vi compatible mode
filetype plugin on " enable plugins
" -- optional --
" auto close options when exiting insert mode
autocmd InsertLeave * if pumvisible() == 0|pclose|endif
set completeopt=menu,menuone
" -- configs --
let OmniCpp_MayCompleteDot = 1 " autocomplete with .
let OmniCpp_MayCompleteArrow = 1 " autocomplete with ->
let OmniCpp_MayCompleteScope = 1 " autocomplete with ::
let OmniCpp_SelectFirstItem = 2 " select first item (but don't insert)
let OmniCpp_NamespaceSearch = 2 " search namespaces in this and included files
let OmniCpp_ShowPrototypeInAbbr = 1 " show function prototype (i.e. parameters) in popup window
" -- ctags --
" map <ctrl>+F12 to generate ctags for current folder:
map <C-F12> :!ctags -R --c++-kinds=+p --fields=+iaS --extra=+q .<CR><CR>
" add current directory's generated tags file to available tags
set tags+=./tags
```

Adding the tags for C++ STL (Standard Template Library)

```
> wget stl-src.zip
> unzip stl-src.zip
> find . -name '*.h' | xargs sed -i 's/__STL_BEGIN_NAMESPACE/namespace std {/'
> find . -name '*.h' | xargs sed -i 's/__STL_END_NAMESPACE/}/'
> ctags -R --c++-kinds=+p --fields=+iaS --extra=+q ./
> mkdir ~/.vim/mytags
> mv tags ~/.vim/mytags/c++stl.tag
> vim ~/.vimrc
vim> set tags+=~/.vim/mytags/c++stl.tag
```

Example code

```
> vim helloworld.cpp
vim> std::[will show autocompletion]
```

Hit **Ctrl + F12** in order to add tags for the current program you're working on

- [Original blog post how-to](http://design.liberta.co.za/articles/code-completion-intellisense-for-cpp-in-vim-with-omnicppcomplete/)

## Tips & tricks

```
# performing edits on multiple files (pipe separates commands)
$ vim -c "argdo %s/ABC/DEF/g | w" *.txt
> vim -c "argdo %s/FOO/BAR/g | update" `grep -l FOO *`
```

Windows users, you might want to check out [Notepad++](http://notepad-plus.sourceforge.net/uk/site.htm) instead of using **vim**.

## External References

- [Vim project homepage](http://www.vim.org/),
- [Cream -- customize vim](http://cream.sourceforge.net/faq.html).
