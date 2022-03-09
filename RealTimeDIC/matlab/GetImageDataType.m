function DataFieldType=GetImageDataType(FileName)
% GetImageDataType - Get image data type.
% 
% AUTHOR:
% 
% Downloaded from Mathowrks
%   Digital Image Correlation and Tracking
%   version 2.1 (7.21 MB) by Melanie Senn
%   https://www.mathworks.com/matlabcentral/fileexchange/50994-digital
%   -image-correlation-and-trackin
% 
%   NOTES:
% Original Function from downloaded code package
% Used to obtain the correct method to convert an image to the
% corresponding graylevel matrix

ImageInfo=imfinfo(FileName);
switch ImageInfo.BitDepth
    case 8
        DataFieldType=@uint8;
    case 16
        DataFieldType=@uint16;
    case 24
        DataFieldType=@uint32;
    case 32
        DataFieldType=@uint32;
    case 64
        DataFieldType=@uint64;
    otherwise
        DataFieldType=[];
end
