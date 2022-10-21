# Dominance-Landscape-Basins

Funkcija `get_dominance_landscape_basins_from_matrix` vrne matrike: 
 - `labeled_pareto_front`:  
 `labeled_pareto_front[x, y] = i`, če je na mestu (x, y) pareto fronta i-tega basin-a, sicer je `labeled_pareto_front[x, y] = 0`
 - `labeled_basins`:  
 `labeled_basins[x, y] = i`, če je mesto (x, y) del i-tega baisina lokalno nedominiranih točk, sicer je `labeled_basins[x, y] = 0`
 - `labeled_slopes`:  
 `labeled_slopes[x, y] = i`, če je labeled_basins[x, y] = i ali od mesta (x, y) vodi pot "navzdol" samo do enega - i-tega baisina, sicer je `labeled_slopes[x, y] = 0`

# Dva primera:

![image](https://user-images.githubusercontent.com/47631585/197243767-438fab66-6685-40ee-8f54-35c82f5f6a34.png)

![image](https://user-images.githubusercontent.com/47631585/197244249-190fc48a-bb93-4350-a53d-84168fb8fb9e.png)
