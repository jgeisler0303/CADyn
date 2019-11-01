files= getnames('..\src\*.cpp');
for i= 1:length(files)
    mex('-c',  '-g', 'CXXFLAGS="$CXXFLAGS -march=native -m64 -std=c++0x"', '-I../src',  '-DEIGEN_DONT_PARALLELIZE', fullfile('../src', files{i}))
end

files= getnames('*.obj');
system(['C:\MATLAB\SupportPackages\R2016a\MW_MinGW_4_9\bin\ar.exe -r "libCADyn.lib" ', sprintf('%s ', files{:})])