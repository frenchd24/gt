 #!/usr/bin/env python
 
def main():

    okList = ['pec','SBc','Sab','SA0','SB(','SA(','SB',\
    'SA',\
    'Sb',\
    'Sc',\
    'S0',\
    'Irr',\
    'dSph',\
    'LINER',\
    'Dwarf',\
    'E5',\
    'Im',\
    'E0',\
    'E+',\
    'BCD',\
    'Sbrst',\
    'Elliptical',\
    'Spiral',\
    'Im_HII',\
    'dwarf']    
      
    MTypes = ['E0 Compact;WR;HII_Sbrst',\
    'x',\
    'S0/a_pec?',\
    "(R'_1)S(r)a",\
    'High_vel._cloud',\
    'DA-star',\
    'High_vel._cloud',\
    'O',\
    'Carbon',\
    'Point_Src_[SDSS]',\
    'Possible_star',\
    'Planetary_nebula']
    
    d_morph = {}
    for MType in MTypes:
        findi = False
        for i in okList:
#             print 'i: ',i
            if MType.find(i) != -1:
                findi = True

            
        if not findi:
            print 'i in MType = ',MType
            print 'MType.find(i) = ', MType.find(i)
                
            if d_morph.has_key(MType):
                d_morph[MType] +=1 
            else:
                d_morph[MType]=1
                
                
#         findi = False
#         for i in okList:
#             if MType.find(i) ==-1:
#                 findi = True
#                 
#         if not findi and MType != 'x':
#             print 'MType = ',MType
#             print 'Name = ',Name  
#             
#             if d_morph.has_key(MType):
#                 d_morph[MType] +=1 
#             else:
#                 d_morph[MType]=1                
                
                
                
                
    print
    print 'd_morph: ',d_morph
    print
                
if __name__ == '__main__':
    main()