function img = normRange(img)
    img = img - min(img(:));
    img = img/max(abs(img(:)));
    return
end