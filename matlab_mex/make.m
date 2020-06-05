function make()
    if (~verLessThan('matlab','9.4'))
        mex -R2018a -output asa_bcp asa_bcp_matlab.cpp
    else
        mex -largeArrayDims -output asa_bcp asa_bcp_matlab.cpp
    end
end