fileName_w='meas_MID00514_FID124501_lb_RosetteSpinEcho_ek_v1a_48x48_te20_w.dat';
kFile='Rosette_traj_48x48.txt';
rose = io_CSIload_twix(fileName_w);
rose.adcDwellTime = 5e-6;
rose = op_CSICombineTime(rose, 'extras');
rose.adcTime = rose.adcDwellTime*[0:(rose.sz(1)-1)];

rose.data=squeeze(rose.data(:,1,:));
rose.sz=size(rose.data);
rose.dims.coils=0;
rose.dims.ky=2;
rose.flags;
rose.flags.coilCombined=1;

%Flatten the k-space data matrix for PSF estimation;
rose.data=ones(rose.sz);

% %Make a gaussian vesrsion:
% x=linspace(-24,24,48);
% y=x;
% [X,Y]=meshgrid(x,y);
% sigma=80;
% G2D=exp(-(((X.^2)+(Y.^2))./(2*sigma)));
% mesh(G2D);
% 
% gaussK=rose;
% gaussK.data=gaussK.data.*G2D;

%Now run the DFT:
PSF_ros = op_CSIFourierTransform(rose, kFile, 'spatial', true, 'spectral', false);
PSF_ros=squeeze(PSF_ros.data(1,:,:));
PSF_ros=PSF_ros/max(abs(PSF_ros(:)));

%Add density compensation:
[dComp_w,sampDensity,kx,ky] = op_CSIPSFCorrection_jn(rose,kFile);

%Now run the DFT on the density compensated trajectory:
PSF_flat = op_CSIFourierTransform(dComp_w, kFile, 'spatial', true, 'spectral', false);
PSF_flat=squeeze(PSF_flat.data(1,:,:));
PSF_flat=PSF_flat/max(abs(PSF_flat(:)));

%Plot the 2D profiles:
figure;mesh(real(PSF_ros));
figure;mesh(real(PSF_flat));

%Make a 1D profile:
%figure;plot([1:48],squeeze(PSF_ros(24,:))/max(PSF_ros(:)),[1:48],squeeze(PSF_flat(24,:))/max(PSF_flat(:)));
figure;plot([1:48],squeeze(PSF_ros(24,:)),[1:48],squeeze(PSF_flat(24,:)));

%Try zeropadding these and making them smooth:
ros_k_cart=fftshift(fft2(PSF_ros));
ros_k_cart_zp=padarray(ros_k_cart,[128,128]);
PSF_ros_zp=ifft2(ifftshift(ros_k_cart_zp));

flat_k_cart=fftshift(fft2(PSF_flat));
flat_k_cart_zp=padarray(flat_k_cart,[128,128]);
PSF_flat_zp=ifft2(ifftshift(flat_k_cart_zp));

%Plot the new 2D Profiles:
figure;mesh(real(PSF_ros_zp));
figure;mesh(real(PSF_flat_zp));

%Make the new 1D profiles:
figure;plot([1:size(PSF_ros_zp,2)],squeeze(PSF_ros_zp(size(PSF_ros_zp,1)/2,:))/max(PSF_ros_zp(:)),[1:size(PSF_flat_zp,2)],squeeze(PSF_flat_zp(size(PSF_ros_zp,1)/2,:))/max(PSF_flat_zp(:)));



%Find the integral under each PSF.  this is roughly equivalent to the SNR
%for each voxel:
SNR_flat=sum(PSF_flat_zp(:));
SNR_rose=sum(PSF_ros_zp(:));

disp(['The SNR of the uncompensated Rosette trajectory is: ' num2str(SNR_rose) '.']);
disp(['The SNR of the flat compensated Rosette trajectory is: ' num2str(SNR_flat) '.']);


