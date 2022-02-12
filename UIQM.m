function uiqm = UIQM(~, c1, c2, c3)
image=imread("C:\Users\piyus\OneDrive\Desktop\Project Phase 2 Update\Output Images\output_7.png");
if ~exist('c1', 'var')
    c1 = 0.0282;
end

if ~exist('c2', 'var')
    c2 = 0.2953;
end

if ~exist('c3', 'var')
    c3 = 3.5753;
end

uicm = UICM(image);
uism = UISM(image);
uiconm = UIConM(image);

uiqm = c1 * uicm + c2 * uism + c3 * uiconm;

end
