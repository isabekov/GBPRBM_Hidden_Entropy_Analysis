function Slider_Elevation_Callback(hObject, ~, ~)
    Elevation = get(hObject, 'Value');
    h_Azimuth = getappdata(hObject, 'h_Slider_Azimuth');
    Azimuth = get(h_Azimuth, 'Value');
    h_Axes_Emp = getappdata(hObject, 'h_Axes_Emp');
    h_Axes_Thr = getappdata(hObject, 'h_Axes_Thr');
    view(h_Axes_Emp, Azimuth, Elevation);
    view(h_Axes_Thr, Azimuth, Elevation);
end