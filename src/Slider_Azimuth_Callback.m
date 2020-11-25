function Slider_Azimuth_Callback(hObject, ~, ~)
    Azimuth = get(hObject, 'Value');
    h_Elevation = getappdata(hObject, 'h_Slider_Elevation');
    Elevation = get(h_Elevation, 'Value');
    h_Axes_Emp = getappdata(hObject, 'h_Axes_Emp');
    h_Axes_Thr = getappdata(hObject, 'h_Axes_Thr');
    view(h_Axes_Emp, Azimuth, Elevation);
    view(h_Axes_Thr, Azimuth, Elevation);
end