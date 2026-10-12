# Flexia

[![Build Status](https://github.com/kutsjuice/Flexia.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/kutsjuice/Flexia.jl/actions/workflows/CI.yml?query=branch%3Amaster)

Библиотека для расчёта и визуализации плоских механизмов: тела, шарниры,
ползуны, заделки, пружины, силы, датчики и приводы.

## Установка

```julia
using Pkg
Pkg.develop(path = "path/to/Flexia")   # или Pkg.add("Flexia")
```

## Визуализация

Визуализация построена на **рецептах Makie**. Пакет зависит только от `Makie`
core, поэтому бэкенд подключается независимо — выберите его в своём скрипте:

```julia
using Flexia
using CairoMakie        # без дисплея: файлы, сервер, CI
# using GLMakie         # интерактивное окно
# using WGLMakie        # браузер/WebGL (запись видео ограничена)
```

Полный цикл — построить сцену, перейти на кадр, записать анимацию:

```julia
fig, ax, p = systemfigure(sys, sol)      # пределы подбираются по траектории
save("frame.png", fig)

record(fig, "out.mp4", 1:5:size(sol, 2); framerate = 40) do t
    p.time = t                           # сцена не пересобирается
end
```

Если нужна своя ось:

```julia
fig = Figure()
ax = Axis(fig[1, 1], aspect = DataAspect())
p = systemplot!(ax, sys, sol; show_forces = false)
```

Интерактивно (GLMakie) — тот же объект графика:

```julia
slider = Slider(fig[1, 2], range = 1:size(sol, 2))
on(slider.value) do t
    p.time = t
end
```

### Визуальные компаньоны

То, **как** выглядит сущность, вынесено из математических типов в отдельные
объекты-компаньоны. Тело — не обязательно «палка»:

```julia
bd = Body2D(1.0, 1.0; length = 0.5)                        # BarVisual, как раньше
bd = Body2D(1.0, 1.0; visual = RectVisual(width = 0.4, height = 0.3))
bd = Body2D(1.0, 1.0; visual = CircleVisual(radius = 0.25))
bd = Body2D(1.0, 1.0; visual = MeshVisual(points = contour))
bd = Body2D(1.0, 1.0; visual = Link2DVisual(length = 1.0, width = 0.25))

set_visualization!(bd, CircleVisual(radius = 0.1))          # поменять на ходу
```

`Link2DVisual` — звено со скруглёнными концами: объединение прямоугольника
`length × width` и двух кружков радиуса `width / 2` с центрами в `±length / 2`.
Центры скруглений совпадают с концами `BarVisual` той же длины, то есть с
типовыми точками крепления шарниров, поэтому такая форма ложится на механизм без
подгонки. Габарит по длине при этом `length + width`.

Готовый пример: `examples/crankshaft_link2d.jl` (копия `crankshaft.jl` с
кривошипом и шатуном в виде `Link2DVisual`).

Компаньон описывает только то, что не меняется во времени: форму, размеры,
прикрепление (смещение и доворот в локальной системе тела) и стиль. Всё, что
зависит от кадра, берётся из решения.

**Стиль.** `color = nothing` (по умолчанию) означает «взять цвет из палитры
Makie» — разные сущности получают разные цвета. Явный цвет слота палитры не
расходует:

```julia
set_visualization!(bd, CircleVisual(radius = 0.1,
                                    style = VisualStyle(color = :red, linewidth = 3)))
```

Цвет может быть функцией от состояния — заложено для раскраски по напряжениям:

```julia
style = VisualStyle(color = (sys, sol, t, entity) -> напряжение > 100 ? :red : :blue)
```

**Не всё должно рисоваться.** У датчиков, приводов и прочих сущностей
компаньона нет — они молча пропускаются. Это легальный случай, а не ошибка:

```julia
visual(sensor) isa NoVisual     # => true, отрисовка не требуется
```

Спрятать отдельный экземпляр, не меняя тип: `VisualStyle(visible = false)`.

### Своя форма тела

Новый тип компаньона плюс два метода — без правок математики и существующих
рецептов:

```julia
struct HexVisual <: AbstractBodyVisual
    radius::Float64
end

function Flexia.body_contour(x, y, θ, v::HexVisual)
    att = Attachment()
    return Point2f[Flexia.to_global(x, y, θ, att, Point2f(v.radius * cos(a), v.radius * sin(a)))
                   for a in range(0, 2π; length = 7)]
end

function Flexia._draw_body!(p::Flexia.BodyPlot, v::HexVisual)
    pts = Makie.lift(p.time) do t
        Flexia.body_contour(p.system[], p.body[], view(p.solution[], :, t))
    end
    Makie.lines!(p, pts)
    return p
end
```

### Порядок отрисовки

`VisualStyle(zorder = ...)`; при равных значениях порядок — как при добавлении:
тела, соединения, силы.

## Структура визуализации

| Файл | Содержимое |
|---|---|
| `src/visuals.jl` | типы компаньонов, `Attachment`, `VisualStyle`, доступ `visual` |
| `src/plotdata.jl` | чистая геометрия по состоянию — без Makie-графиков, тестируется без бэкенда |
| `src/recipes.jl` | рецепты Makie; зависимость только от `Makie`, не от бэкенда |

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
