To install ADMB, you need to install a C++ compiler as well as ADMB. I suggest that you install gcc.The following are
instructions for installation with Windows. For other operating systems see (http://www.admb-project.org/downloads) and 
(http://www.admb-project.org/tools/gcc/). 

Windows Instructions:

Install gcc 4.5 from: http://www.admb-project.org/tools/gcc/gcc452-win32.zip/view
Put in c:/MinGW

Install ADMB 11: http://admb-project.googlecode.com/files/admb-11-mingw-gcc4.5-32bit.exe
Put in C:/admb

I use the following function in R to setup R2admb to access ADMB rather than adding to my path so gcc versions
with Rtools don't conflict. 

prepare_admb=function()
{
    Sys.setenv(PATH = paste("c:/admb/bin;c:admb/utilities;c:/MinGW/bin;", 
        Sys.getenv("PATH"), sep = ";"))
    Sys.setenv(ADMB_HOME = "c:/admb")
    invisible()
}

To use different locations you'll need to change the values used above


