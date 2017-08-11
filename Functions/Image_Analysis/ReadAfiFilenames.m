function fileNames=ReadAfiFilenames(afiFileName)
    
    xmldata=xmlread(afiFileName);
    
    import javax.xml.xpath.*
    factory = XPathFactory.newInstance;
    xpath = factory.newXPath;
    
    % compile and evaluate the XPath Expression
    expression = xpath.compile('ImageList/Image/Path');
    nodeList= expression.evaluate(xmldata, XPathConstants.NODESET);
    
    fileNames=cell(nodeList.getLength,1);
    for i = 1:nodeList.getLength
        node = nodeList.item(i-1);
        fileNames{i}=char(node.getFirstChild.getNodeValue);
    end
end