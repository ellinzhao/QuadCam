% Multi-sensor DiffuserCam via 
 
%   Add NET assembly if it does not exist
%   May need to change specific location of library
 
%{
asm = System.AppDomain.CurrentDomain.GetAssemblies;
if ~any(arrayfun(@(n) strncmpi(char(asm.Get(n-1).FullName), ...
        'uEyeDotNet', length('uEyeDotNet')), 1:asm.Length))
    NET.addAssembly(...
        'C:\Program Files\IDS\uEye\Develop\DotNet\signed\uEyeDotNet.dll');
end
%   Create camera object handle
cam = uEye.Camera;
%   Open 1st available camera
%   Returns if unsuccessful
if ~strcmp(char(cam.Init), 'SUCCESS')
    error('Could not initialize camera');
end
%   Set display mode to bitmap (DiB)
if ~strcmp(char(cam.Display.Mode.Set(uEye.Defines.DisplayMode.DiB)), ...
        'SUCCESS')
    error('Could not set display mode');
end
%   Set colormode to 8-bit RAW
if ~strcmp(char(cam.PixelFormat.Set(uEye.Defines.ColorMode.SensorRaw8)), ...
        'SUCCESS')
    error('Could not set pixel format');
end
%   Set trigger mode to software (single image acquisition)
if ~strcmp(char(cam.Trigger.Set(uEye.Defines.TriggerMode.Software)), 'SUCCESS')
    error('Could not set trigger format');
end
%   Allocate image memory
[ErrChk, img.ID] = cam.Memory.Allocate(true);
if ~strcmp(char(ErrChk), 'SUCCESS')
    error('Could not allocate memory');
end
%   Obtain image information
[ErrChk, img.Width, img.Height, img.Bits, img.Pitch] ...
    = cam.Memory.Inquire(img.ID);
if ~strcmp(char(ErrChk), 'SUCCESS')
    error('Could not get image information');
end
%   Acquire image
if ~strcmp(char(cam.Acquisition.Freeze(true)), 'SUCCESS')
    error('Could not acquire image');
end
%   Extract image
[ErrChk, tmp] = cam.Memory.CopyToArray(img.ID); 
if ~strcmp(char(ErrChk), 'SUCCESS')
    error('Could not obtain image data');
end
%   Reshape image
img.Data = reshape(uint8(tmp), [img.Width, img.Height, img.Bits/8]);
%   Draw image
himg = imshow(img.Data, 'Border', 'tight');
%   Acquire & draw 100 times
for n=1:100
    %   Acquire image
    if ~strcmp(char(cam.Acquisition.Freeze(true)), 'SUCCESS')
        error('Could not acquire image');
    end
      %   Extract image
      [ErrChk, tmp] = cam.Memory.CopyToArray(img.ID); 
      if ~strcmp(char(ErrChk), 'SUCCESS')
          error('Could not obtain image data');
      end
      %   Reshape image
      img.Data = reshape(uint8(tmp), [img.Width, img.Height, img.Bits/8]);
      %   Draw image
      set(himg, 'CData', img.Data);
      drawnow;
  end
%   Close camera
if ~strcmp(char(cam.Exit), 'SUCCESS')
    error('Could not close camera');
end
 
%}

% Image save folder
folder = 'T:\Deshler\XYStage_Diffuser_Profiling\Diffuser_profile_imgs_RPI\';
mkdir(folder);

% User should orient XY stage at top left of sensor.
imaq; 
origin = [mmc.getXPosition('XYStage'), mmc.getYPosition('XYStage')];

sensor_dim = [7.2 5.4];
w = 7.2; % sensor width (mm)
h = 5.4; % sensor height (mm)
s = 200; % sensor spacing (mm) 
diffuser_dim = [300 300]; 
offset = .5 * (diffuser_dim - 2*sensor_dim  - s);

% Setup camera object
cam = camera();
cam.Brightness = 50;
cam.Contrast = 0;
cam.Saturation = -100;
cam.Sharpness = 0;
cam.ExposureMode = 'night';
cam.ExposureCompensation = 0;
cam.AWBMode = 'off';

for n = 0:4
    i = mod(n, 2);
    j = n / 2;
    x = i * (sensor_dim(1) + s) + offset(1);
    y = j * (sensor_dim(2) + s) + offset(2);
    mmc.setXYPosition(x, y);
    mmc.waitForDevice('XYStage');
    capture(cam, folder, n);
end
 

function capture(cam, folder, i)    
    img = snapshot(cam);
    imwrite(img, [folder,'img' ,int2str(i),'.png']);
    clear img;    
end
