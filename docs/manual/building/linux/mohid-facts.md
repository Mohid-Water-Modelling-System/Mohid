# Mohid facts

Linux is great to quickly get some trivia related with the MOHID software suite. These numbers contemplate only the MOHIDWater, MOHIDLand, MOHIDBase1 and MOHIDBase2 code and omit the suite of MOHID related tools built in Fortran or Visual Basic.

## Code related facts (as of 2010-04-21)

How many **keywords** are there in MohidWater (including Base1 and Base2 libraries)? **3436**

```
> for i in `find | grep MOHID | grep -v Land | grep -v River | grep F90`; do cat $i | grep -re "call GetData" ; done | wc -l
```

---

And in MohidWater and MohidLand together? **4082**

```
> for i in `find | grep MOHID | grep F90`; do cat $i | grep -re "call GetData" ; done | wc -l
```

---

How many **lines of code** is there in MohidWater and MohidLand (including Base1 and Base2)? **703820**

```
> for i in `find | grep MOHID | grep F90`; do cat $i; done | wc -l
```

---

And without empty lines? **476330**

```
> for i in `find | grep MOHID | grep F90`; do cat $i | grep "\w"; done | wc -l
```

---

And without comment-lines? **400500**

```
> for i in `find | grep MOHID | grep F90`; do cat $i | grep "\w" | grep -v "^ *\!"; done | wc -l
```

---

How many parallel **openmp** zones is there in the code? **1428**

```
> for i in `find | grep MOHID | grep F90`; do cat $i | grep "OMP PARALLEL"; done | wc -l
```

---

How many subroutines is there in the code? **6466**

```
> for i in `find | grep MOHID | grep F90`; do cat $i | grep -r "^ *subroutine"; done | wc -l
```

---

How many **ifs** are there in the code? **40366**

```
> for i in `find | grep MOHID | grep F90`; do cat $i | grep -r "^ *if ("; done | wc -l
```

---

How many **calls** are there in the code? **29770**

```
> for i in `find | grep MOHID | grep F90`; do cat $i | grep -r "^ *call "; done | wc -l
```

---

How many **types** are there in the code? **4994**

```
> for i in `find | grep MOHID | grep F90`; do cat $i | grep -r "^ *type "; done | wc -l
```

---

How many **allocates** are there in the code? **1860**

```
> for i in `find | grep MOHID | grep F90`; do cat $i | grep -r "^ *allocate "; done | wc -l
```

---

And **deallocates**? **1322**

```
> for i in `find | grep MOHID | grep F90`; do cat $i | grep -r "^ *deallocate "; done | wc -l
```

Could there be memory leaks? Maybe not, as one deallocate command may free several allocated arrays... It would be a good thing enforce one deallocate command per allocate command. That would make a quick check for memory leaks in the program...

---

How many **do loops** are there in the code? **8096**

```
> for i in `find | grep MOHID | grep F90`; do cat $i | grep -r "^ *do "; done | wc -l
```

---

## Related links

- [Codeplex](../../archive/web/codeplex.md)
- [Regexp](../../archive/linux/regexp.md)
- [Bash](../../archive/linux/bash.md)
- [grep](../../archive/linux/grep.md)
