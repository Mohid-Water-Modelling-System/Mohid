# Sam2p

**sam2p** is a little utility that allows to convert between image file formats. Better than [convert](../../building/windows/convert.md) (generated eps filesizes are one-hundredfold smaller).

## Example

Converting from jpg to eps:

```
> sam2p file.jpg EPS: file.eps
```

Batch convert from jpg to eps (windows):

```
> for %i in (*.jpg) do sam2p %i EPS: %~ni.eps
```

## See also

- [Convert](../../building/windows/convert.md)

## External links

- [sam2p google code homepage.](http://code.google.com/p/sam2p/)
