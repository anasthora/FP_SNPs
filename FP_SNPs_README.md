# FP_SNPs

Краткое описание: скрипт для выявления референсного и альтернативного аллелей в файле формата «#CHROM<TAB>POS<TAB>ID<TAB>allele1<TAB>allele2»

## Содержание
- [Подготовка](#подготовка)
- [Использование](#использование)
- [Лицензия](#лицензия)

## Подготовка
Получение предобработанного файла FP_SNPs_10k_GB38_twoAllelsFormat.tsv на основе исходного FP_SNPs.txt.

запуск команды в терминале:
```bash
awk 'BEGIN{FS=OFS="\t"} NR==1{print "#CHROM","POS","ID","allele1","allele2"; next} $2==23{next} {print "chr"$2, $4, "rs"$1, $5, $6}' FP_SNPs.txt > FP_SNPs_10k_GB38_twoAllelsFormat.tsv
```
В файле FP_SNPs_10k_GB38_twoAllelsFormat.tsv:
- отсутствуют координаты по GRCh37 (GB37_position)
- переименованы столбцы и изменен их порядок;
- добавлены префиксы chr и rs;
- удалены варианты с X-хромосомы.

## Использование
Скрипт resolve_ref_alt.py на вход принимает:
- список маркеров FP_SNPs_10k_GB38_twoAllelsFormat.tsv с аллелями в формате "allele1","allele2"
- референсный геном GRCh38.d1.vd1.fa 

В результате работы скрипт выдает:
- список маркеров с указанием референсного и альтернативного аллеля
- список маркеров, для которых не найдено совпадение с референсом (с указанием ошибки)
- лог-файл

запуск крипта в терминале:
```bash
python3 resolve_ref_alt.py \
    -i FP_SNPs_10k_GB38_twoAllelsFormat.tsv \
    -f GRCh38.d1.vd1.fa \
    -o FP_SNPs_10k_GB38_final.tsv  \
    --err FP_SNPs_10k_GB38_errors.tsv \
    --log FP_SNPs_10k_GB38.log
```

## Лицензия
без лицензии

```bash
git clone https://github.com/anasthora/FP_SNPs.git
cd FP_SNPs
