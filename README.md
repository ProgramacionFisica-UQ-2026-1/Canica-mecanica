# README — Simulador de Canica sobre Superficie Topográfica

## Descripción general

Este proyecto implementa un simulador físico en Python que modela el movimiento de una **canica (partícula puntual)** sobre diferentes superficies topográficas bidimensionales. La partícula se mueve bajo la influencia de:

- **Gravedad** (acelera según la pendiente del terreno)
- **Rozamiento viscoso** (lineal con la velocidad)
- **Topografía variable** (el gradiente de la superficie determina la dirección de la fuerza)

El sistema resuelve numéricamente las ecuaciones diferenciales del movimiento mediante el método `RK45` de `SciPy`, calcula las energías cinética, potencial y total, y genera visualizaciones gráficas de la trayectoria y la evolución energética. Además, incluye un **análisis avanzado** que ajusta parámetros físicos (rozamiento y gravedad) a partir de datos sintéticos con ruido, utilizando `curve_fit`.

El programa está diseñado para ser **interactivo por consola**, guiando al usuario paso a paso: selección de terreno, posición inicial, previsualización del mapa, posibilidad de cambiar de terreno o posición, y finalmente ingreso de velocidad, rozamiento y gravedad para ejecutar la simulación.

---

## Características principales

- ✅ Simulación dinámica de partículas sobre superficies arbitrarias definidas por funciones matemáticas.
- ✅ Implementación **orientada a objetos** con:
  - Clase abstracta `SistemaFisico` (define la interfaz).
  - Clase `Terreno` (almacena la superficie, calcula gradiente y energía potencial).
  - Clase `Particula` (hereda de `SistemaFisico` e implementa la integración).
- ✅ **Eventos automáticos** durante la integración con `solve_ivp`:
  - Detiene la simulación cuando la velocidad cae por debajo de un umbral (la canica se para).
  - Detiene la simulación si la partícula sale del dominio `[-7, 7]` (evita cálculos fuera de rango).
- ✅ Cálculo y visualización de **energías cinética, potencial y total**.
- ✅ **Gráficas profesionales**:
  - Mapa topográfico (contourf) con la trayectoria superpuesta (punto verde = inicio, azul = fin).
  - Evolución temporal de las energías en una misma figura.
- ✅ **Previsualización interactiva del terreno** antes de la simulación final:
  - Tras elegir el terreno y la posición inicial, se muestra un mapa con un punto rojo indicando la posición.
  - El usuario debe **presionar Enter** para cerrar la ventana del gráfico y continuar.
  - A continuación, aparece un menú para:
    1. Continuar con la simulación (usar ese terreno y posición).
    2. Cambiar de terreno (elegir otro mapa).
    3. Cambiar solo la posición (mantener el mismo terreno).
  - Así se puede explorar diferentes configuraciones sin reiniciar el programa.
- ✅ **Validación robusta de entradas**: rangos, números, etc.
- ✅ **Análisis avanzado con `curve_fit`**:
  - Genera datos sintéticos (trayectoria real + ruido).
  - Ajusta los parámetros de rozamiento y gravedad para que la simulación se aproxime a los datos ruidosos.
  - Calcula MSE y MAE, mostrando la incertidumbre de los parámetros.

---

## Estructura del proyecto

```text
proyecto/
├── simulador.py          # Código fuente completo
└── README.md             # Este archivo
