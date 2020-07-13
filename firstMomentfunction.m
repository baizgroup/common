function firstMoment = firstMomentfunction(NORM_poly_FTIR_Bkg_Corr, interpAxis, input1, input2)

 for n=1
 %Defines the start and end wavenumbers of the peak intensities
peakWindow(:,n) = NORM_poly_FTIR_Bkg_Corr(input1:input2,n);  
FreqWindow(:,n) = interpAxis(input1 : input2);

%total area of the peak
totalArea(:,n) =sum(peakWindow(:,n));
%normalizes each intensity value by the total area 
areaNorm = peakWindow(:,n)./totalArea(:,n) ;

%creates a set of intensity and wavenumber
normSpec = [ FreqWindow(:,n) , areaNorm ];
%increases scale of intensity by multiplying it by freq


freqAreaProd(:,n) = peakWindow(:,n).* FreqWindow(:,n);
tops(:,n) = sum(freqAreaProd(:,n));
bottoms(:,n) = sum(peakWindow(:,n));
centerArea(:,n) = tops(:,n)./bottoms(:,n)
disp(num2str(centerArea(:,n), '%.2f'))
firstMoment(:,n) = centerArea(:,n)
 end
end 