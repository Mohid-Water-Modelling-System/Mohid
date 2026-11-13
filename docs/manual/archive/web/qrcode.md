# Qrcode

> **Archived page.** This content is kept for historical reference only. It may be out of date, incomplete or unrelated to the current MOHID documentation, and it is not maintained.

[Qrcode](http://en.wikipedia.org/wiki/QR_Code) is a non-proprietary encoding of information into 2d codes. 2d codes are bidimensional matrices that encode information much like the traditional barcodes you find in supermarket products. Widely spread in Japan for some years now, it allows users with cameras on their mobile phones to scan the 2d code and extract the information. If it's a hyperlink, then the scanner software may open the phone's web browser and open the hyperlink page.

## Hello world example

[http://farm4.static.flickr.com/3093/2824480513\_cb08e32b3b\_m.jpg](http://farm4.static.flickr.com/3093/2824480513_cb08e32b3b_m.jpg)

### Python

You can easily use the [python](../linux/python.md) library [PyQrCodec](http://www.pedemonte.eu/wiki/index.php/Documentation:PyQrCodec) (depends on the [PIL](http://www.pythonware.com/products/pil/) library)

```
> sudo yum install opencv opencv-devel opencv-python
```

```
> wget http://effbot.org/downloads/Imaging-1.1.6.tar.gz
> tar -xvf Imaging-1.1.6.tar.gz
> cd Imaging-1.1.6
> python setup.py build
> sudo python setup.py install
```

```
> wget http://www.pedemonte.eu:81/pyqr/files/PyQrcodec_Linux.tar.gz
> tar -xvf PyQrcodec_Linux.tar.gz
> cd PyQrcodec
> python setup.py build
> sudo python setup.py install
```

`[python] >>> import PyQrCodec >>> size, image = PyQrCodec.encode("Hello world!") >>> image.save("helloworld.jpg")`

## External links

- [QR code in wikipedia](http://en.wikipedia.org/wiki/QR_Code)

### Encoders/decoders

- [Semapedia](http://semapedia.org)
- [Kaywa](http://qrcode.kaywa.com/)
- [Neoreader](http://www.neoreader.com/code.html)
