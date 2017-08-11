function [vertices,annoText]=Get_Vertices_From_XML(xml_filename,nDim)
    if(nargin<2)
       nDim=2; 
    end
    xmlres=parseXML(xml_filename);
    
    annotation_layers=find(strcmp({xmlres.Children.Name},'Annotation'));
    number_of_annotations=length(annotation_layers);
    
    vertices=cell(number_of_annotations,1);
    annoText=cell(number_of_annotations,1);
    
    for annotation_number=1:number_of_annotations
        
        annotation=xmlres.Children(annotation_layers(annotation_number));
        region_list_layer=find(strcmp({annotation.Children.Name},'Regions'),1);
        region_list=annotation.Children(region_list_layer);
        region_list_layers=find(strcmp({region_list.Children.Name},'Region'));
        number_of_regions=length(region_list_layers);
        vertices{annotation_number}=cell(number_of_regions,1);
        annoText{annotation_number}=cell(number_of_regions,1);
        
        for region_number=1:number_of_regions
            region=region_list.Children(region_list_layers(region_number));
            textLayer=find(strcmp({region.Attributes.Name},'Text'));
            annoText{annotation_number}{region_number}=region.Attributes(textLayer).Value;
            vertices_layer=find(strcmp({region.Children.Name},'Vertices'),1);
            vertice_list=region.Children(vertices_layer);
            vertex_layers=find(strcmp({vertice_list.Children.Name},'Vertex'));
            number_of_vertices=length(vertex_layers);
            vertmat=zeros(nDim,number_of_vertices);
            temp=[vertice_list.Children(vertex_layers).Attributes];
            vertmat(:)=str2double({temp(:).Value});
            vertices{annotation_number}{region_number}=vertmat';
        end
        
    end
end
