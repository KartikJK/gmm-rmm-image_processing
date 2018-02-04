load feature
load trainfea
load target
c = knnclassify(feature, trainfea, target);
if c==1
   
  msgbox('Calcified Plaque')
end
if c==2
 msgbox('Fibrotic Plaque')
end
if c==3
    
     msgbox('Lipidic Plaque')

end