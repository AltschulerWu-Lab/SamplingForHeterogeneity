function tumor_outline=Load_Tumor_Outline(outlineXmlFile,image_type)
if(nargin<2)
    image_type='IF';
end

if(strcmp(image_type,'IF'))
    vertices=Get_Vertices_From_XML(outlineXmlFile);
    xvals=vertices{1}{1}(:,1);
    yvals=vertices{1}{1}(:,2);
    xvals(end)=xvals(1);
    yvals(end)=yvals(1);
    tumor_outline=[xvals,yvals];
else
    vertices=load(obj.params.tumor_outline_filename);
    tumor_outline=vertices.tissue_vertices{1}';
    tumor_outline=bsxfun(@minus,tumor_outline,min(tumor_outline));
end

end