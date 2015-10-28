@echo off
call c:\MinGW32_4.6.1\mingwvars.bat

del rate_sse3.obj
del rate_sse.obj
gcc rate_sse3.c  -c -O2 -march=core2 -mno-ssse3 -mfpmath=sse -mstackrealign -ffast-math -D WIN32 -D NDEBUG -D _WINDOWS -o rate_sse3.obj
gcc rate_sse.c   -c -O2 -march=core2 -msse -mno-sse2 -mno-sse3 -mno-ssse3 -mstackrealign -ffast-math -D WIN32 -D NDEBUG -D _WINDOWS -o rate_sse.obj
copy _divdi3.461 _divdi3.obj

cd ffmpeg
del fft_mmx.obj
del fft_sse.obj
yasm -f win32 --prefix=_ -o fft_mmx.obj fft_mmx.asm
gcc fft_sse.c -c -O2 -Wall -mfpmath=sse -mstackrealign -ffast-math -march=core2 -msse -mno-sse2 -mno-sse3 -mno-ssse3 -D WIN32 -D NDEBUG -D _WINDOWS -o fft_sse.obj
cd ..

cd fft4g
del fft4g_sse3.obj
gcc fft4g_sse3.c -c -O2 -march=core2 -mno-ssse3 -mfpmath=sse -mstackrealign -ffast-math -D WIN32 -D NDEBUG -D _WINDOWS -o fft4g_sse3.obj
cd ..

pause
