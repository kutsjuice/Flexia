# Руководство пользователя Flexia

Исходники документации на Quarto. Собирается в сайт (HTML), DOCX и PDF.

## Сборка

```powershell
quarto render .            # все форматы → _book/
quarto preview .           # предпросмотр сайта
```

Результат: `_book/index.html`, `_book/Flexia-manual.pdf`,
`_book/Flexia-manual.docx`.

Для PDF нужен `xelatex` (MiKTeX/TeX Live). Шрифты с кириллицей заданы в
`_quarto.yml` (`Times New Roman`, `Arial`, `Consolas`); на не-Windows системах
замените их на `DejaVu Serif`, `DejaVu Sans`, `DejaVu Sans Mono`.

## Публикация

```powershell
.\deploy-pages.ps1
```

Скрипт собирает документацию и публикует `_book/` в ветку `pages` remote
`origin`. Ветка `pages` — сгенерированная. GitHub Pages должен быть настроен
на источник «ветка `pages`, каталог `/`».

Подробности — в главе `publishing.qmd`.
