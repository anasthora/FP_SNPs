#!/usr/bin/env python3
"""
Скрипт для выявления референсного и альтернативного аллелей во входном файле формата «#CHROM<TAB>POS<TAB>ID<TAB>allele1<TAB>allele2».
Выходные файлы:
  1) основной файл *_final.tsv
  2) файл ошибок *_errors.tsv - строки, где ни allele1, ни allele2 не совпали с REF из FASTA
  3) лог-файл со статистикой 

Перед запуском выявления референсного и альтернативного аллелей проверяет наличие индекса .fai.
"""

import argparse
import logging
import os
import sys
import pysam


def setup_logging(log_path: str, verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    fmt = "%(asctime)s [%(levelname)s] %(message)s"
    datefmt = "%Y-%m-%d %H:%M:%S"

    # Логи в stderr
    stream_handler = logging.StreamHandler(sys.stderr)
    stream_handler.setFormatter(logging.Formatter(fmt, datefmt=datefmt))

    # Логи в файл
    file_handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    file_handler.setFormatter(logging.Formatter(fmt, datefmt=datefmt))

    logging.basicConfig(level=level, handlers=[stream_handler, file_handler])


def parse_args():
    parser = argparse.ArgumentParser(
        prog="resolve_ref_alt.py",
        description=(
            "По файлу '#CHROM\\tPOS\\tID\\tallele1\\tallele2' и референсному "
            "геному определяет REF и ALT. "
            "Пишет основной файл, файл ошибок и лог-файл."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Пример:\n"
            "  %(prog)s -i variants.txt -f GRCh38.d1.vd1.fa \\\n"
            "      -o out.tsv --err err.tsv --log run.log\n"
        ),
    )
    parser.add_argument("-i", "--input", required=True, metavar="FILE",
                        help="Входной файл со столбцами #CHROM POS ID allele1 allele2")
    parser.add_argument("-f", "--fasta", required=True, metavar="FILE",
                        help="Референсный FASTA (должен иметь .fai индекс)")
    parser.add_argument("-o", "--output", required=True, metavar="FILE",
                        help="Основной выходной файл")
    parser.add_argument("--err", required=True, metavar="FILE",
                        help="Файл ошибок (строки без совпадения с REF)")
    parser.add_argument("--log", required=True, metavar="FILE",
                        help="Лог-файл со статистикой")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Подробный вывод (DEBUG)")
    return parser.parse_args()


def ensure_fasta_index(fasta_path: str) -> None:
    """Проверяет наличие .fai; при отсутствии — создаёт через pysam.faidx."""
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA не найден: {fasta_path}")
    fai_path = fasta_path + ".fai"
    if not os.path.exists(fai_path):
        logging.warning("Индекс .fai не найден, создаю: %s", fai_path)
        pysam.faidx(fasta_path)
        if not os.path.exists(fai_path):
            raise IOError(f"Не удалось создать индекс: {fai_path}")
    logging.info("Индекс FASTA найден: %s", fai_path)


def check_header(line: str, line_no: int) -> None:
    expected = ["#CHROM", "POS", "ID", "allele1", "allele2"]
    fields = line.rstrip("\r\n").split("\t")
    if fields[: len(expected)] != expected:
        raise ValueError(
            f"Строка {line_no}: неверный заголовок.\n"
            f"  Ожидается: {expected}\n"
            f"  Получено:  {fields[: len(expected)]}"
        )
    logging.info("Заголовок проверен: %s", "\t".join(fields[: len(expected)]))


def resolve_ref_alt(fasta: pysam.FastaFile, chrom: str, pos: int,
                    a1: str, a2: str):
    """
    Возвращает (ref, alt) или None, если ни один аллель не совпал с геномом.
    """
    if chrom not in fasta.references:
        raise ValueError(f"Хромосома '{chrom}' отсутствует в FASTA")
    if pos < 1 or pos > fasta.get_reference_length(chrom):
        raise ValueError(
            f"Позиция {pos} вне диапазона хромосомы {chrom} "
            f"(длина {fasta.get_reference_length(chrom)})"
        )

    ref_base = fasta.fetch(chrom, pos - 1, pos).upper()
    a1u, a2u = a1.upper(), a2.upper()

    m1 = (a1u == ref_base)
    m2 = (a2u == ref_base)

    if m1 and not m2:
        return a1, a2
    if m2 and not m1:
        return a2, a1
    if m1 and m2:
        # Гомозигота: оба совпадают с референсом
        return a1, a2
    return None  # ни один не совпал


def main():
    args = parse_args()
    setup_logging(args.log, args.verbose)
    logging.info("Запуск: %s", " ".join(sys.argv))

    try:
        ensure_fasta_index(args.fasta)
    except (FileNotFoundError, IOError) as e:
        logging.error("Ошибка индекса FASTA: %s", e)
        sys.exit(2)

    if not os.path.exists(args.input):
        logging.error("Входной файл не найден: %s", args.input)
        sys.exit(2)

    n_ok = 0
    n_err = 0
    n_skipped = 0

    try:
        with pysam.FastaFile(args.fasta) as fasta, \
             open(args.input, "r", encoding="utf-8", newline="") as fin, \
             open(args.output, "w", encoding="utf-8", newline="\n") as fout, \
             open(args.err, "w", encoding="utf-8", newline="\n") as ferr:

            # Заголовок
            header_line = fin.readline()
            if not header_line:
                raise ValueError("Входной файл пуст.")
            check_header(header_line, 1)
            fout.write("#CHROM\tPOS\tID\tREF\tALT\n")
            # В файле ошибок сохраним исходные поля + причину
            ferr.write("#CHROM\tPOS\tID\tallele1\tallele2\treason\n")

            for line_no, raw in enumerate(fin, start=2):
                line = raw.rstrip("\r\n")
                if not line:
                    continue

                fields = line.split("\t")
                if len(fields) < 5:
                    logging.warning("Строка %d: ожидалось 5 колонок, получено %d — пропуск.",
                                    line_no, len(fields))
                    n_skipped += 1
                    continue

                chrom, pos_s, vid, a1, a2 = fields[:5]

                # Парсинг позиции
                try:
                    pos = int(pos_s)
                except ValueError:
                    logging.warning("Строка %d: нечисловая POS=%r — пропуск.",
                                    line_no, pos_s)
                    n_skipped += 1
                    continue

                # Разрешение REF/ALT
                try:
                    result = resolve_ref_alt(fasta, chrom, pos, a1, a2)
                except ValueError as e:
                    logging.warning("Строка %d: %s — строка в файл ошибок.",
                                    line_no, e)
                    ferr.write(f"{chrom}\t{pos_s}\t{vid}\t{a1}\t{a2}\t{e}\n")
                    n_err += 1
                    continue

                if result is None:
                    ref_base = fasta.fetch(chrom, pos - 1, pos).upper()
                    reason = (f"ни allele1={a1!r}, ни allele2={a2!r} "
                              f"не совпали с REF={ref_base!r}. Возможно, комплементарная цепь.")
                    logging.warning("Строка %d: %s", line_no, reason)
                    ferr.write(
                        f"{chrom}\t{pos_s}\t{vid}\t{a1}\t{a2}\t{reason}\n"
                    )
                    n_err += 1
                    continue

                ref, alt = result
                fout.write(f"{chrom}\t{pos_s}\t{vid}\t{ref}\t{alt}\n")
                n_ok += 1

                if (n_ok + n_err) % 10000 == 0:
                    logging.info("Обработано: ok=%d, err=%d", n_ok, n_err)

    except Exception as e:
        logging.error("Критическая ошибка: %s", e)
        sys.exit(1)

    # Итоговая статистика
    logging.info("Готово.")
    logging.info("Записано в основной файл (%s): %d", args.output, n_ok)
    logging.info("Записано в файл ошибок (%s): %d", args.err, n_err)
    logging.info("Пропущено некорректных строк: %d", n_skipped)

    print(f"OK={n_ok}\tERR={n_err}\tSKIPPED={n_skipped}")


if __name__ == "__main__":
    main()