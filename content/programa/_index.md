---
linkTitle: Programa
title: "Programa del curso"
type: book
date: "`r as.Date(file.mtime(knitr::current_input()))`"
toc: true
weight: 10
---

{{< toc hide_on="xl" >}}

## Objetivos

* Conocer la teoría sobre la que se fundamentan los métodos para la estimación de relaciones empíricas y la inferencia usando datos de sección cruzada y de panel.
* Diseñar estrategias econométricas usando los modelos adecuados de acuerdo con la pregunta de investigación.
* Emplear software para estimar los modelos econométricos apropiados de acuerdo con la naturaleza de los datos disponibles.
* Conocer e implementar buenas prácticas en el uso de software, orientadas a la transparencia y la replicabilidad en la investigación.
* Conocer los métodos que se emplean en la investigación económica actual.

## Referencias 

El curso se basa en los siguientes textos:

* (MHE) Angrist, J.D. y Pischke, J.S. (2013). Mostly *Harmless Econometrics: An Empiricists Companion*. Princeton University Press.
* \* (CT) Cameron, A.C. y Trivedi, P.K. (2005). *Microeconometrics: Methods and applications*. Oxford University Press.
*	Hayashi, F. (2000). *Econometrics*. Princeton University Press.
*	\* (W) Wooldridge, J.M. (2010). *Econometric analysis of cross section and panel data*. Segunda edición, MIT Press.

## Contenido temático

Unidad 1. Introducción
+ Revisión de fundamentos de estadística y regresión lineal
+ Teoría asintótica
+ Estimadores extremos
+ Máxima verosimilitud
+ Mínimos cuadrados no lineales
+ Inferencia estadística


Unidad 2. Modelos de variable dependiente no continua  
+ Modelos de variable dependiente binaria  
+ Modelos de conteo  

Unidad 3. Modelos de selección
+ Modelo de Tobit
+ Modelo de Heckman

Unidad 4. Endogeneidad
+ Variables instrumentales
+ Método generalizado de momentos
+ Estimación con instrumentos débiles

Unidad 5. Datos de panel
+ Modelos de efectos fijos y de efectos aleatorios
+ Estimadores between y within
+ Estimadores de primeras diferencias y de efectos aleatorios
+ Errores estándar agrupados

Unidad 6. Extensiones
+ Bootstrap
+ Modelos de riesgo y de sobrevivencia
+ Regresión por cuantiles
+ Modelos de panel con endogeneidad
+ Modelos de panel no lineales
+ Métodos semi paramétricos y no paramétricos

## Evaluación del curso

+ Examen parcial: 30%.
+ Examen final acumulativo: 45%
+ Tareas (4): 20% (5% cada una)
+ Exposición: 5%

## Tareas

Cuatro tareas teórico-prácticas. Las tareas deben entregarse de manera individual, pero se recomienda ampliamente colaborar en grupos de estudio. Las tareas deberán entregarse en Teams antes de la fecha y hora señalada. No se aceptarán tareas fuera de tiempo. Por favor, no comprima los archivos en carpetas comprimidas. Las tareas deberán contener dos archivos:
Un primer documento de respuestas donde se incluyan las respuestas a las preguntas teóricas y conceptuales. Este documento debe estar en formato pdf y debe ser generado usando un software de procesamiento de textos científicos, por ejemplo, usando los leguajes LaTeX o Markdown. En este documento también se deben incluir las respuestas a preguntas sobre conclusiones que se desprenden de las secciones prácticas. Por ejemplo, si una pregunta pide obtener la media de la variable x en cierta base de datos, entonces el documento de respuestas debe incluir la pregunta y respuesta correspondiente: “la media de la variable x es 32.6”. En este documento también deberán incluirse las tablas y gráficas que se soliciten.

Un segundo archivo deberá contener el código replicable usado para generar los resultados de la sección práctica. El código debe también crear las tablas y gráficas solicitadas. Los archivos de código se verificarán para comprobar su replicabilidad de manera aleatoria.

## Software

R será el paquete standard usado en las sesiones prácticas. Más aún, el uso de cualquier software es aceptado siempre que se cumplan con los requisitos de replicabilidad y reportes de las tareas y exámenes. Habrá una sesión especial para introducir el uso de RMarkdown para la generación de reportes científicos.

## Exámenes
+ Examen parcial: jueves 6 de octubre de 2022 en el horario de clase.
+ Examen final: martes 13 de diciembre de 2022 en horario por definir.

## Exposición

Cada alumno realizará una presentación de uno de los artículos aplicados marcados con negritas en la lista de lecturas. Cada presentación deberá ser de máximo 15 minutos y deberá incluir el contenido que el presentador considere más relevante. La presentación deberá abordar, mínimamente: 1) el problema a investigar, 2) la metodología empleada, 3) la relación entre la metodología y la teoría vista en el curso, 4) los datos empleados, 5) los principales resultados, y 6) una crítica sobre la validez y las conclusiones del estudio.
