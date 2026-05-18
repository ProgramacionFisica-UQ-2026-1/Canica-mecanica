# README — Simulador de Canica sobre Superficie Topográfica

## Descripción

Este proyecto implementa un simulador físico en Python para modelar el movimiento de una canica sobre diferentes superficies topográficas bidimensionales bajo la influencia de:

* Gravedad
* Rozamiento viscoso
* Pendiente del terreno

El sistema resuelve numéricamente las ecuaciones diferenciales del movimiento utilizando el método `RK45` de `SciPy`, calcula energías mecánicas y genera visualizaciones gráficas de la trayectoria y evolución energética.

El proyecto también incluye un análisis avanzado mediante ajuste de parámetros usando `curve_fit`.

---

## Características principales

* Simulación dinámica de partículas sobre superficies arbitrarias.
* Implementación orientada a objetos.
* Uso de clases abstractas (`ABC`).
* Integración numérica de EDOs con `solve_ivp`.
* Detección automática de eventos:

  * Detención por velocidad mínima.
  * Salida del dominio.
* Cálculo de:

  * Energía cinética
  * Energía potencial
  * Energía mecánica total
* Visualización gráfica:

  * Trayectoria sobre mapa topográfico
  * Evolución temporal de energías
* Ajuste de parámetros físicos mediante `curve_fit`.
* Validación robusta de entradas del usuario.

---

## Estructura del proyecto

```text
proyecto/
│
├── simulador.py
├── README.md
```

---

## Requisitos

Instalar las siguientes librerías:

```bash
pip install numpy matplotlib scipy
```

---

## Librerías utilizadas

| Librería          | Uso                                    |
| ----------------- | -------------------------------------- |
| `numpy`           | Operaciones numéricas                  |
| `matplotlib`      | Gráficas y visualización               |
| `scipy.integrate` | Resolución de ecuaciones diferenciales |
| `scipy.optimize`  | Ajuste de parámetros                   |
| `abc`             | Clases abstractas                      |

---

## Modelo físico

La aceleración de la partícula se calcula a partir del gradiente del terreno y el rozamiento:

a_x=-g\frac{\partial z}{\partial x}-\mu v_x

a_y=-g\frac{\partial z}{\partial y}-\mu v_y

donde:

* ( g ): gravedad
* ( \mu ): coeficiente de rozamiento
* ( z(x,y) ): altura del terreno

La energía cinética se calcula mediante:

E_c=\frac{1}{2}m(v_x^2+v_y^2)

La energía potencial gravitatoria:

E_p=mgz(x,y)

---

## Terrenos disponibles

El simulador incluye cinco superficies topográficas predefinidas:

1. Montaña clásica
2. Plano inclinado con montículos
3. Dos picos y valle central
4. Silla de montar
5. Montaña rugosa

Cada terreno se define mediante funciones matemáticas dependientes de (x) y (y).

---

## Arquitectura del código

### 1. Clase abstracta `SistemaFisico`

Define la interfaz base del sistema físico:

* `derivada()`
* `simular()`

---

### 2. Clase `Terreno`

Representa la superficie topográfica.

Funciones principales:

* cálculo de altura
* gradiente numérico
* energía potencial

---

### 3. Clase `Particula`

Hereda de `SistemaFisico`.

Responsable de:

* resolver ecuaciones diferenciales
* calcular aceleraciones
* ejecutar simulaciones físicas

---

### 4. Eventos de simulación

#### `PararPorVelocidad`

Finaliza la simulación cuando:

```text
|v| < umbral
```

#### `LimiteDominio`

Detiene la simulación si la partícula sale del rango:

```text
[-7, 7]
```

---

## Ejecución

Ejecutar el archivo principal:

```bash
python simulador.py
```

---

## Flujo de uso

El programa solicita:

1. Selección del terreno.
2. Posición inicial.
3. Velocidad inicial.
4. Rozamiento.
5. Gravedad.

Posteriormente:

* ejecuta la simulación,
* calcula estadísticas,
* genera gráficas,
* realiza ajuste de parámetros.

---

## Resultados generados

### Estadísticas

El programa muestra:

* tiempo total
* posición inicial y final
* velocidad inicial y final
* rapidez media
* aceleración media
* distancia recorrida
* pérdida de energía

---

### Gráficas

#### Trayectoria topográfica

Muestra:

* superficie del terreno
* recorrido de la canica
* punto inicial
* punto final

#### Energías

Grafica:

* energía cinética
* energía potencial
* energía total

---

## Ajuste de parámetros (`curve_fit`)

El proyecto incluye un módulo avanzado que:

1. Genera datos simulados con ruido.
2. Ajusta:

   * rozamiento
   * gravedad
3. Evalúa:

   * MSE
   * MAE

Esto permite estudiar técnicas básicas de inferencia y optimización numérica.

---

## Conceptos aplicados

Este proyecto integra conceptos de:

* Mecánica clásica
* Métodos numéricos
* Ecuaciones diferenciales ordinarias
* Optimización
* Programación orientada a objetos
* Análisis energético
* Simulación computacional

---

## Posibles mejoras

* Animación en tiempo real.
* Exportación de resultados.
* Interfaz gráfica.
* Terrenos cargados desde archivos externos.
* Simulación 3D.
* Integración con motores físicos.
* Paralelización de simulaciones.


