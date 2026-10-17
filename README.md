# Esther

A Next-Generation CXX 20 Hardware Identification Kit based on floating-point computations via CRT.

## Exclamation-Mark

This project is a POC. It is not very stable yet, but has been tested successfully on different machines.
This only supports the MSVC compiler and is specifically made for Windows.
For MSVC, the needed compiler flags are strict

```cl /std:c++20 /GS- program.cpp /link /NODEFAULTLIB /ENTRY:entry /SUBSYSTEM:CONSOLE kernel32.lib```

Every program also needs this:

```extern "C" int _fltused = 0;```

in its source files, because the floating-point code references _fltused in a non-CRT link.

## Question-mark

The ```harware_identification()``` combines two hashes into a 64-bit value. The first hashes the Library's arithmetic over 4096 fixed generated samples, and that part is the same on any IEEE 754 machine. The second hashes ```exp```, ```log```, ```sin``` and ```pow``` from ```ucrtbase.dll``` which is loaded at runtime. Since it changes on different CPU builds and Windows versions, a unique fingerprint of the hardware is computed.

## Add-mark

To use this library just include the header ```esther.hpp```
and call the corresponding function ```esther::hardware_identification();``` it returns a unsigned 64 bit integer

## Tilde-mark 

Expanding on this, you can use this to 

``` Lock the fingerprint via a licensing system to a client and let it compute a part of the decryption key on the fly based on the library to decrypt incoming networking packets from the server. By doing this, you're locking a key to this specific device. If someone were to share the loader now, the decryption would fail. ```

This still doesn't fully replace the traditional path of hardware identification (for now, who knows, maybe in the future). It's smarter to layer this on top of or combine it with the traditional print.

## Add

- Claude has assisted in the research
- Inspiration is taken from browser fingerprinting libs, so I have made it available for the desktop too.
