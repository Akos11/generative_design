# generative_design

## Beüzemelés

A környezet teljesen a 3D geometria tárgyban kapott keretre épül. A beüzemelése azzal megegyezik.

## Teszt-modellek

Látható, hogy a teszthez használt modellek a `testmodels/` mappában vannak. Egy adott teszt-modell több különböző fájlból épül fel. Ezért készítettem egy összefogó .txt -t hogy könnyebb legyen a betöltés.
A betöltéshez az elindított tesztprogramon belül a `File/Open` menün belül a load-in.txt kell kiválsztani, amely betölti az összes modellt megfelelően. `(testmodels/2D_01/load-in.txt)`

## Használat

A jelenlegi kódban az eredeti környezetben használt gyorsbillentyűk működnek. pl. ` W - wireframed `
Az új funkciókhoz tartozó billentyűk:
- 1 - a kényszer alakzatok megjelenítése/eltűntetése
- 2 - a háromszögek szegmentálása az alapján, hogy melyek tartoznak egyes kényszer alakzatokhoz