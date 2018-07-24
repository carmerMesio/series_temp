En este documento se explica el material adjunto proporcionado.


* Output.pdf 
Trabajo en formato pdf

* Output.Rmd
Markdown del trabajo

* animate_with_exo.html
Performance del algoritmo en el conjunto test considerando var. exógenas
Se adjunta código para replicar (animate.R)
Es necesario tener la carpeta animate_images, js y css en el mismo directorio que 
el .html

* Archivo csv

Archivo de datos original Quarter Hourly Electrical Energy_01-01-2017_31-12-2017.csv

Scripts
-------

* chosing.R
Función implementada para evaluar los pesos óptimos para el algoritmo ponderado AR, ets.
Parámetros: ch= {0,1}  0 si no se quiere var. exo, 1 al contrario

* animate.R
Implementación que permite ver las animación para el mejor resultado en el conjunto test.

* scrap.R

Implementación que permite obtener los datos festivos del año 2017.
Fuente: https://www.calendario-365.es/dias-festivos/2017.html

* Files Rdata
Se adjuntan 3 archivos que incluyen un dataframe con la base de datos original
, variables exógenas y estimación de la función Fafitstats.



