# Flexia

[![Build Status](https://github.com/kutsjuice/Flexia.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/kutsjuice/Flexia.jl/actions/workflows/CI.yml?query=branch%3Amaster)

Библиотека для расчёта и визуализации плоских механизмов: тела, шарниры,
ползуны, заделки, пружины, силы, датчики и приводы.

## Установка

```julia
using Pkg
Pkg.develop(path = "path/to/Flexia")   # или Pkg.add("Flexia")
```

## Документация

Полное руководство пользователя — установка, модель, связи, силы, приводы,
датчики, расчёт и визуализация:

- **HTML:** <https://kutsjuice.github.io/Flexia/>
- **PDF и DOCX:** кнопки скачивания на той же странице.

Исходники руководства лежат в [`docs/`](docs/) и собираются
[Quarto](https://quarto.org) в HTML, DOCX и PDF:

```powershell
quarto render docs                 # все форматы
.\docs\deploy-pages.ps1            # собрать и опубликовать в ветку pages
```

Визуализация описана в главе
[«Визуализация»](https://kutsjuice.github.io/Flexia/visualization.html)
руководства.

## Тесты

```julia
using Pkg; Pkg.test()
```

Тесты геометрии и компаньонов работают без установленного бэкенда; рендер
проверяется headless на `CairoMakie`.

## Примеры

В каталоге `examples/`:

```bash
julia --project=examples examples/pendulum.jl
```
