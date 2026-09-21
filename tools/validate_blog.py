#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Read-only validator for the blog post contract.

The validator never writes or normalizes files. It reports contract gaps so
existing posts can be migrated deliberately, while --strict can be used for
new or ready-to-publish posts.
"""

from __future__ import print_function

import argparse
import datetime as dt
import html as html_std
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple
from urllib.parse import unquote, urlparse

try:
    import yaml
except ImportError:  # pragma: no cover - environment guard
    sys.stderr.write("需要 PyYAML，请先安装：python -m pip install pyyaml\n")
    sys.exit(2)


ROOT = Path(__file__).resolve().parents[1]
TAXONOMY_PATH = ROOT / "docs" / "blog-taxonomy.yml"
AUTHORS_PATH = ROOT / "_data" / "authors.yml"
POST_DIR = ROOT / "_posts"
DRAFT_DIR = ROOT / "_drafts"

POST_NAME_RE = re.compile(
    r"^(\d{4}-\d{2}-\d{2})-([\u4e00-\u9fffA-Za-z0-9-]+)\.md$"
)
DRAFT_NAME_RE = re.compile(r"^([\u4e00-\u9fffA-Za-z0-9-]+)\.md$")
DATE_RE = re.compile(
    r"^(\d{4})-(\d{2})-(\d{2})(?:\s+(\d{1,2}):(\d{2}):(\d{2})\s+([+-]\d{4}))?$"
)
MARKDOWN_IMAGE_START_RE = re.compile(r"!\[[^\]]*\]\(")
HTML_IMAGE_RE = re.compile(
    r"<img\b[^>]*\bsrc\s*=\s*([\"'])(.*?)\1", re.IGNORECASE | re.DOTALL
)


@dataclass(frozen=True)
class Finding:
    level: str
    code: str
    path: Path
    message: str


class Report:
    def __init__(self) -> None:
        self.findings: List[Finding] = []

    def add(self, level: str, code: str, path: Path, message: str) -> None:
        finding = Finding(level, code, path, message)
        if finding not in self.findings:
            self.findings.append(finding)

    def error_count(self) -> int:
        return sum(1 for item in self.findings if item.level == "ERROR")

    def warning_count(self) -> int:
        return sum(1 for item in self.findings if item.level == "WARN")

    def print(self, checked_count: int, mode: str) -> None:
        ordered = sorted(
            self.findings,
            key=lambda item: (display_path(item.path), item.level, item.code),
        )
        for item in ordered:
            print(
                "[{level}] {path} [{code}] {message}".format(
                    level=item.level,
                    path=display_path(item.path),
                    code=item.code,
                    message=item.message,
                )
            )
        print("")
        print(
            "检查 {count} 篇文章；模式：{mode}；错误 {errors}；警告 {warnings}".format(
                count=checked_count,
                mode=mode,
                errors=self.error_count(),
                warnings=self.warning_count(),
            )
        )


def display_path(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT)).replace("\\", "/")
    except ValueError:
        return str(path)


def load_yaml_file(path: Path) -> object:
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def load_taxonomy() -> Dict[str, object]:
    if not TAXONOMY_PATH.exists():
        raise RuntimeError("找不到分类文件：{}".format(display_path(TAXONOMY_PATH)))

    data = load_yaml_file(TAXONOMY_PATH)
    if not isinstance(data, dict):
        raise RuntimeError("分类文件必须使用 YAML 映射作为顶层结构")

    categories = data.get("categories")
    if not isinstance(categories, dict) or not categories:
        raise RuntimeError("分类文件缺少 categories")

    for top, children in categories.items():
        if not isinstance(top, str) or not isinstance(children, list):
            raise RuntimeError("categories 必须是“一级分类: [二级分类]”结构")
        if not children or not all(isinstance(item, str) and item for item in children):
            raise RuntimeError("一级分类 {} 的二级分类不能为空".format(top))

    tags = data.get("tags")
    if not isinstance(tags, list) or not tags:
        raise RuntimeError("分类文件缺少 tags")
    if len(tags) != len(set(tags)):
        raise RuntimeError("tags 白名单存在重复项")

    aliases = data.get("tag_aliases", {})
    if not isinstance(aliases, dict):
        raise RuntimeError("tag_aliases 必须是映射")

    return data


def load_authors() -> Set[str]:
    if not AUTHORS_PATH.exists():
        raise RuntimeError("找不到作者文件：{}".format(display_path(AUTHORS_PATH)))
    data = load_yaml_file(AUTHORS_PATH)
    if not isinstance(data, dict):
        raise RuntimeError("_data/authors.yml 必须使用 YAML 映射")
    return {str(key) for key in data.keys()}


def discover_files(paths: Sequence[str]) -> List[Path]:
    if paths:
        found: List[Path] = []
        for raw_path in paths:
            candidate = Path(raw_path)
            if not candidate.is_absolute():
                candidate = Path.cwd() / candidate
            candidate = candidate.resolve()
            if candidate.is_dir():
                found.extend(sorted(candidate.rglob("*.md")))
            elif candidate.suffix.lower() == ".md":
                found.append(candidate)
        return sorted(set(found), key=lambda item: str(item).lower())

    found = []
    for directory in (POST_DIR, DRAFT_DIR):
        if directory.exists():
            found.extend(sorted(directory.rglob("*.md")))
    return sorted(set(found), key=lambda item: str(item).lower())


def infer_kind(path: Path) -> str:
    parts = set(path.parts)
    if "_posts" in parts:
        return "post"
    if "_drafts" in parts:
        return "draft"
    if POST_NAME_RE.match(path.name):
        return "post"
    return "draft"


def read_front_matter(path: Path) -> Tuple[Optional[Dict[str, object]], str, Optional[str]]:
    text = path.read_text(encoding="utf-8")
    if text.startswith("\ufeff"):
        text = text[1:]

    lines = text.splitlines()
    if not lines or lines[0].strip() != "---":
        return None, text, "缺少 front matter 起始标记"

    end_index = None
    for index in range(1, len(lines)):
        if lines[index].strip() == "---":
            end_index = index
            break
    if end_index is None:
        return None, text, "front matter 缺少结束标记"

    front_matter_text = "\n".join(lines[1:end_index])
    body = "\n".join(lines[end_index + 1 :])
    try:
        front_matter = yaml.safe_load(front_matter_text)
    except yaml.YAMLError as exc:
        return None, body, "front matter YAML 无法解析：{}".format(exc)

    if not isinstance(front_matter, dict):
        return None, body, "front matter 必须是 YAML 映射"
    return front_matter, body, None


def normalize_date(value: object) -> Tuple[Optional[dt.datetime], Optional[str]]:
    if isinstance(value, dt.datetime):
        return value, value.strftime("%Y-%m-%d %H:%M:%S %z")
    if isinstance(value, dt.date):
        return None, value.isoformat()

    text = str(value).strip()
    match = DATE_RE.match(text)
    if not match:
        return None, text

    timezone_text = match.group(7)
    timezone = None
    if timezone_text:
        offset_hours = int(timezone_text[1:3])
        offset_minutes = int(timezone_text[3:5])
        offset = dt.timedelta(hours=offset_hours, minutes=offset_minutes)
        if timezone_text[0] == "-":
            offset = -offset
        timezone = dt.timezone(offset)

    try:
        parsed = dt.datetime(
            int(match.group(1)),
            int(match.group(2)),
            int(match.group(3)),
            int(match.group(4) or 0),
            int(match.group(5) or 0),
            int(match.group(6) or 0),
            tzinfo=timezone,
        )
    except ValueError:
        return None, text
    return parsed, text


def count_han(text: str) -> int:
    return len(re.findall(r"[\u4e00-\u9fff]", text))


def count_english_words(text: str) -> int:
    return len(re.findall(r"[A-Za-z]+(?:'[A-Za-z]+)?", text))


def strip_fenced_code(text: str) -> str:
    output: List[str] = []
    in_fence = False
    for line in text.splitlines():
        if re.match(r"^\s*(```|~~~)", line):
            in_fence = not in_fence
            output.append("")
            continue
        output.append("" if in_fence else line)
    return "\n".join(output)


def extract_image_sources(text: str) -> List[str]:
    content = strip_fenced_code(text)
    sources: List[str] = []
    for match in MARKDOWN_IMAGE_START_RE.finditer(content):
        start = match.end()
        depth = 1
        index = start
        while index < len(content) and depth:
            character = content[index]
            if character == "(":
                depth += 1
            elif character == ")":
                depth -= 1
            index += 1
        if depth == 0:
            sources.append(content[start : index - 1])
    sources.extend(match.group(2) for match in HTML_IMAGE_RE.finditer(content))
    return sources


def resolve_local_image(source: str, document_path: Path) -> Optional[Path]:
    value = html_std.unescape(source).strip()
    if value.startswith("<") and value.endswith(">"):
        value = value[1:-1].strip()
    title_match = re.match(r"^(.*?)(?:\s+[\"'][^\"']*[\"'])$", value)
    if title_match:
        value = title_match.group(1).strip()

    parsed = urlparse(value)
    if parsed.scheme.lower() in {"http", "https", "data", "mailto", "tel"}:
        return None
    if value.startswith("#"):
        return None

    clean_path = unquote(parsed.path)
    if clean_path.startswith("/"):
        return (ROOT / clean_path.lstrip("/")).resolve()
    return (document_path.parent / clean_path).resolve()


def is_under(path: Path, directory: Path) -> bool:
    try:
        path.relative_to(directory)
        return True
    except ValueError:
        return False


def category_alias_map(taxonomy: Dict[str, object]) -> Dict[Tuple[str, ...], Tuple[str, ...]]:
    result: Dict[Tuple[str, ...], Tuple[str, ...]] = {}
    aliases = taxonomy.get("category_aliases", [])
    if not isinstance(aliases, list):
        return result
    for item in aliases:
        if not isinstance(item, dict):
            continue
        source = item.get("source")
        target = item.get("target")
        if isinstance(source, list) and isinstance(target, list):
            result[tuple(str(value) for value in source)] = tuple(
                str(value) for value in target
            )
    return result


def check_filename(
    path: Path,
    kind: str,
    front_matter: Optional[Dict[str, object]],
    report: Report,
) -> None:
    if kind == "post":
        match = POST_NAME_RE.match(path.name)
        if not match:
            report.add(
                "ERROR",
                "filename",
                path,
                "发布文章文件名必须为 YYYY-MM-DD-中文语义标题.md，且不能包含空格、括号或全角标点",
            )
            return

        filename_date = match.group(1)
        if front_matter is None or "date" not in front_matter:
            return
        parsed_date, _ = normalize_date(front_matter.get("date"))
        if parsed_date is not None and parsed_date.strftime("%Y-%m-%d") != filename_date:
            report.add(
                "ERROR",
                "date_filename_mismatch",
                path,
                "front matter 日期 {0} 与文件名日期 {1} 不一致".format(
                    parsed_date.strftime("%Y-%m-%d"), filename_date
                ),
            )
        return

    if not DRAFT_NAME_RE.match(path.name):
        report.add(
            "ERROR",
            "draft_filename",
            path,
            "草稿文件名只能使用中文、ASCII 字母数字和连字符，不能带日期",
        )


def check_front_matter(
    path: Path,
    kind: str,
    front_matter: Dict[str, object],
    taxonomy: Dict[str, object],
    authors: Set[str],
    report: Report,
) -> Optional[bool]:
    required = [
        "title",
        "description",
        "author",
        "categories",
        "tags",
        "pin",
    ]
    if kind == "post":
        required.insert(3, "date")

    for field in required:
        if field not in front_matter or front_matter.get(field) in (None, ""):
            report.add("ERROR", "required_field", path, "缺少必填字段：{}".format(field))

    title = front_matter.get("title")
    if title is not None and (not isinstance(title, str) or not title.strip()):
        report.add("ERROR", "title", path, "title 必须是非空字符串")

    description = front_matter.get("description")
    if description is not None:
        if not isinstance(description, str) or not description.strip():
            report.add("ERROR", "description", path, "description 必须是非空字符串")
        elif "\n" in description:
            report.add("ERROR", "description", path, "description 不能使用多行文本")
        else:
            han_count = count_han(description)
            if han_count:
                if han_count < 50 or han_count > 120:
                    report.add(
                        "WARN",
                        "description_length",
                        path,
                        "中文 description 当前 {0} 个汉字，建议控制在 50–120 个".format(
                            han_count
                        ),
                    )
            else:
                words = count_english_words(description)
                if words < 15 or words > 80:
                    report.add(
                        "WARN",
                        "description_length",
                        path,
                        "英文 description 当前 {0} 个单词，建议控制在 15–80 个".format(
                            words
                        ),
                    )

    author = front_matter.get("author")
    if author is not None:
        if not isinstance(author, str) or not author.strip():
            report.add("ERROR", "author", path, "author 必须是非空字符串")
        elif author not in authors:
            report.add(
                "ERROR",
                "unknown_author",
                path,
                "作者 {} 未在 _data/authors.yml 中登记".format(author),
            )

    if "date" in front_matter and front_matter.get("date") not in (None, ""):
        parsed_date, date_text = normalize_date(front_matter.get("date"))
        if parsed_date is None:
            if kind == "post":
                report.add(
                    "ERROR",
                    "date",
                    path,
                    "date 必须为 YYYY-MM-DD HH:MM:SS +0800，当前为 {!r}".format(
                        date_text
                    ),
                )
        elif parsed_date.tzinfo is None:
            if kind == "post":
                report.add("ERROR", "date_timezone", path, "date 缺少 +0800 时区")
        elif parsed_date.strftime("%z") != "+0800":
            report.add(
                "ERROR",
                "date_timezone",
                path,
                "date 必须使用 +0800，当前为 {}".format(parsed_date.strftime("%z")),
            )

    categories = front_matter.get("categories")
    categories_ok = False
    subcategory: Optional[str] = None
    category_tree = taxonomy["categories"]
    assert isinstance(category_tree, dict)

    if not isinstance(categories, list):
        report.add("ERROR", "categories", path, "categories 必须是两级 YAML 数组")
    elif len(categories) != 2:
        report.add(
            "ERROR",
            "categories",
            path,
            "categories 必须恰好包含一个一级分类和一个二级分类",
        )
    elif not all(isinstance(item, str) and item.strip() for item in categories):
        report.add("ERROR", "categories", path, "categories 的每一项都必须是非空字符串")
    else:
        top, sub = categories[0], categories[1]
        subcategory = sub
        children = category_tree.get(top)
        if not isinstance(children, list):
            report.add("ERROR", "category_top", path, "未知的一级分类：{}".format(top))
        elif sub not in children:
            report.add(
                "ERROR",
                "category_sub",
                path,
                "{} 下没有二级分类：{}".format(top, sub),
            )
        else:
            categories_ok = True

        alias_map = category_alias_map(taxonomy)
        source = tuple(str(item) for item in categories)
        if source in alias_map and source != alias_map[source]:
            report.add(
                "WARN",
                "category_migration",
                path,
                "旧分类可迁移为：{}".format(
                    ", ".join(alias_map[source])
                ),
            )

    tags = front_matter.get("tags")
    tag_whitelist = set(str(item) for item in taxonomy["tags"])
    aliases = taxonomy.get("tag_aliases", {})
    assert isinstance(aliases, dict)

    if not isinstance(tags, list):
        report.add("ERROR", "tags", path, "tags 必须是 YAML 数组")
    else:
        if len(tags) < 2 or len(tags) > 6:
            report.add(
                "ERROR",
                "tag_count",
                path,
                "tags 每篇必须为 2–6 个，当前为 {} 个".format(len(tags)),
            )

        seen_tags: Set[str] = set()
        for tag in tags:
            if not isinstance(tag, str) or not tag.strip():
                report.add("ERROR", "tag_value", path, "标签必须是非空字符串")
                continue
            if "，" in tag:
                report.add(
                    "ERROR",
                    "tag_punctuation",
                    path,
                    "标签 {!r} 含中文逗号，请拆成多个标签并使用 ASCII 逗号".format(
                        tag
                    ),
                )
            if tag in seen_tags:
                report.add("ERROR", "tag_duplicate", path, "标签重复：{}".format(tag))
            seen_tags.add(tag)

            if tag in aliases:
                report.add(
                    "ERROR",
                    "tag_alias",
                    path,
                    "标签 {!r} 应替换为规范写法：{}".format(tag, aliases[tag]),
                )
            elif tag not in tag_whitelist:
                report.add(
                    "ERROR",
                    "tag_unknown",
                    path,
                    "标签未登记在 taxonomy 白名单：{}".format(tag),
                )

            if categories_ok and subcategory is not None and tag == subcategory:
                report.add(
                    "ERROR",
                    "tag_equals_category",
                    path,
                    "标签不能与当前二级分类相同：{}".format(tag),
                )

    pin = front_matter.get("pin")
    if not isinstance(pin, bool):
        report.add("ERROR", "pin", path, "pin 必须是 true 或 false")
    return pin if isinstance(pin, bool) else None


def check_images(path: Path, body: str, report: Report) -> None:
    for source in extract_image_sources(body):
        resolved = resolve_local_image(source, path)
        if resolved is None:
            continue
        if not resolved.exists():
            report.add(
                "ERROR",
                "image_missing",
                path,
                "找不到图片：{}".format(source),
            )
            continue

        images_root = ROOT / "assets" / "images"
        if is_under(resolved, images_root):
            relative_parts = resolved.relative_to(images_root).parts
            if relative_parts:
                actual_directory = relative_parts[0]
                expected_directory = path.stem
                if actual_directory != expected_directory:
                    report.add(
                        "ERROR",
                        "image_directory",
                        path,
                        "图片目录为 assets/images/{0}/，应与文章目录 assets/images/{1}/ 一致".format(
                            actual_directory, expected_directory
                        ),
                    )


def check_code_blocks(path: Path, body: str, report: Report) -> None:
    fence_lines = [
        line
        for line in body.splitlines()
        if re.match(r"^\s*(```|~~~)", line)
    ]
    if len(fence_lines) % 2 != 0:
        report.add("ERROR", "code_fence", path, "代码围栏数量为奇数，可能没有闭合")

    without_code = strip_fenced_code(body)
    if without_code.count("$$") % 2 != 0:
        report.add("ERROR", "math_delimiter", path, "$$ 数学公式分隔符数量为奇数")


def check_math_switches(
    path: Path,
    front_matter: Dict[str, object],
    body: str,
    report: Report,
) -> None:
    math = front_matter.get("math", False)
    mermaid = front_matter.get("mermaid", False)

    if not isinstance(math, bool):
        report.add("ERROR", "math_switch", path, "math 必须是 true 或 false")
    if not isinstance(mermaid, bool):
        report.add("ERROR", "mermaid_switch", path, "mermaid 必须是 true 或 false")

    without_code = strip_fenced_code(body)
    has_math = bool(
        re.search(r"\$\$|(?<!\$)\$[^$\n]+\$(?!\$)", without_code)
    )
    has_mermaid = bool(re.search(r"(?m)^\s*```mermaid\b", body))

    if isinstance(math, bool):
        if not math and has_math:
            report.add(
                "ERROR",
                "math_switch",
                path,
                "正文包含数学公式，但 math 不是 true",
            )
        elif math and not has_math:
            report.add("WARN", "math_switch", path, "math: true，但正文没有检测到公式")

    if isinstance(mermaid, bool):
        if not mermaid and has_mermaid:
            report.add(
                "ERROR",
                "mermaid_switch",
                path,
                "正文包含 Mermaid 图，但 mermaid 不是 true",
            )
        elif mermaid and not has_mermaid:
            report.add(
                "WARN",
                "mermaid_switch",
                path,
                "mermaid: true，但正文没有检测到 Mermaid 图",
            )


def check_orphan_asset_directories(
    files: Sequence[Path],
    taxonomy: Dict[str, object],
    report: Report,
) -> None:
    images_root = ROOT / "assets" / "images"
    if not images_root.exists():
        return

    expected = {path.stem for path in files}
    exempt = {
        str(item)
        for item in taxonomy.get("asset_exempt_directories", [])
        if isinstance(item, str)
    }

    for child in sorted(images_root.iterdir(), key=lambda item: item.name):
        if not child.is_dir() or child.name in exempt or child.name in expected:
            continue
        report.add(
            "ERROR",
            "orphan_asset_dir",
            child,
            "图片目录没有对应文章，可能是孤立目录或命名不一致",
        )


def check_file(
    path: Path,
    kind: str,
    taxonomy: Dict[str, object],
    authors: Set[str],
    report: Report,
) -> Optional[bool]:
    front_matter, body, front_matter_error = read_front_matter(path)
    if front_matter_error:
        report.add("ERROR", "front_matter", path, front_matter_error)
        front_matter = None

    pin_value: Optional[bool] = None
    if front_matter is not None:
        pin_value = check_front_matter(
            path, kind, front_matter, taxonomy, authors, report
        )
        check_math_switches(path, front_matter, body, report)

    check_filename(path, kind, front_matter, report)
    check_images(path, body, report)
    check_code_blocks(path, body, report)
    return pin_value


def main() -> int:
    parser = argparse.ArgumentParser(
        description="只读校验博客 front matter、分类、标签和资源引用。"
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--audit",
        action="store_true",
        help="审计现有文章；发现错误时仍返回成功，适合迁移前盘点。",
    )
    mode.add_argument(
        "--strict",
        action="store_true",
        help="严格校验新文章或待发布文章；发现错误时返回失败。",
    )
    parser.add_argument(
        "paths",
        nargs="*",
        help="可选的文件或目录；默认检查 _posts 和 _drafts。",
    )
    args = parser.parse_args()

    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

    try:
        taxonomy = load_taxonomy()
        authors = load_authors()
    except (OSError, RuntimeError, yaml.YAMLError) as exc:
        print("配置错误：{}".format(exc))
        return 2

    files = discover_files(args.paths)
    report = Report()
    pin_count = 0

    for path in files:
        if path.name == ".placeholder":
            continue
        kind = infer_kind(path)
        pin_value = check_file(path, kind, taxonomy, authors, report)
        if pin_value is True:
            pin_count += 1

    if not args.paths:
        check_orphan_asset_directories(files, taxonomy, report)

    pin_limit = int(taxonomy.get("pin_limit", 2))
    if pin_count > pin_limit:
        report.add(
            "ERROR",
            "pin_limit",
            ROOT,
            "置顶文章数量为 {}，超过全站上限 {}".format(pin_count, pin_limit),
        )

    checked_count = len([path for path in files if path.name != ".placeholder"])
    mode_name = "strict" if args.strict else "audit"
    report.print(checked_count, mode_name)

    if args.strict and report.error_count():
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
