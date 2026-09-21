# 博客写作约定（Chirpy Starter / Jekyll）

站点：`https://tangb.rrt01.cn`，主题 `jekyll-theme-chirpy`，`lang: zh`。
推送 `main` 后由 `.github/workflows/pages-deploy.yml` 构建并部署到 GitHub Pages。
CI 中带 `htmlproofer`，**内部链接或图片失效会让构建失败**，所以改完文章必须校验。

## 助手在本仓库的默认动作

1. 写新文章：先落到 `_drafts/`，用《Front matter 模板》；用户确认后再移到 `_posts/` 并加 `YYYY-MM-DD-` 前缀。
2. 分类和标签只使用《分类与标签规范》中列出的值。要新增一级/二级分类，先问用户。
3. 每次改动后运行 `pwsh -File tools/check-posts.ps1` 校验文件名、front matter、日期、分类标签、图片链接。
4. 不要手改 `_site/`（构建产物，已 gitignore）、`assets/lib/`（主题自带）、`assets/js/dist/`。
5. 改动前先看 `git status`；用户的工作区改动一律保留。

## 文件命名

- 文章：`_posts/YYYY-MM-DD-标题.md`，文件名日期必须等于 front matter 的 `date` 日期。
  **没有 `YYYY-MM-DD-` 前缀的 md 文件 Jekyll 完全不构建，线上看不到**（本仓库曾踩过这个坑）。
- 标题中不要出现空格、`()`、`：`、`/` 等符号，因为它们会进入 permalink `/posts/:title/`。
  历史上已存在的 `2025-01-03-call variants.md`、`2025-06-24-docker 简单教程.md` 保留不动，避免旧链接失效；新文章一律不带空格。
- 图片：`assets/images/<与文章同名（含日期前缀）的目录>/xxx.png`，正文引用写
  `![描述](../assets/images/<目录>/xxx.png)`。
  **图片目录名里不要有括号**，Markdown 链接解析会在 `)` 处截断（2025-05-01-pyMC 曾因此踩坑）。
- 未完成的文章放 `_drafts/`，文件名不必带日期；预览需要 `bundle exec jekyll s --drafts`。

## Front matter 模板

```yaml
---
title: 文章标题
description: 一句话摘要，会用于 SEO 和首页卡片，建议 40-90 字。
author: tangb
date: 2026-01-01 09:33:00 +0800
categories: [一级分类, 二级分类]
tags: [标签1, 标签2, 标签3]
pin: false
math: true
mermaid: true
---
```

- `math: true` 只在正文有 `$...$` 公式时开；`mermaid: true` 只在有 mermaid 图时开。
- `pin: true` 会把文章钉在首页顶部，目前旧文章几乎全部是 `true`。**新文章默认 `pin: false`**，只在用户明确要求时置顶。
- 协作/转载他人内容时把 `author` 写成原作者名，不要默认 `tangb`。

## 分类与标签规范

分类回答「属于哪个方向」，标签回答「涉及哪些具体概念」，两者不要重复。

### 允许的一级 / 二级分类

| 一级分类 | 二级分类 | 覆盖内容 |
| --- | --- | --- |
| 生物信息学 | 测序与变异 | 测序原理、比对、变异检测、质控 |
| 生物信息学 | 基因组分析 | 基因预测、注释、非编码 RNA、重复序列、比较基因组 |
| 统计方法 | 回归模型 | 线性回归、加权、广义线性模型 |
| 统计方法 | 蒙特卡洛方法 | MCMC、抽样、贝叶斯计算、pyMC |
| 统计遗传学 | 孟德尔随机化 | MR 原理、IVW、敏感性分析 |
| 统计遗传学 | 遗传率估计 | LDSC、遗传率、遗传相关 |
| 软件工具 | 环境管理 | Docker、conda、依赖打包 |
| 学术写作 | 论文各节 | Abstract / Introduction / Methods / Results / Discussion |

`categories` 必须是上表中的两级，不要出现三级（例：`[软件工具, 环境管理, docker]` 是错的）。

### 标签规则

- 每篇 3-5 个，用英文半角逗号分隔，**禁止全角逗号 `，`**（会被解析成一个假标签，例：`基因组分析，汇总`）。
- 用具体术语：`MCMC`、`GATK`、`LDSC`、`IVW`、`异方差`；不要用「科普」「汇总」「笔记」这类无检索价值的词。
- 不要重复分类名（`categories: [生物信息学, 测序与变异]` 时，标签里不要再出现「生物信息学」）。
- 同一个概念固定一种写法（统一用 `MCMC`、`GWAS`、`孟德尔随机化`，不要中英混用多种拼写）。
- 学术写作类每篇固定带 `论文结构` + 本节英文名 + 该节专属标签。

## 发布前检查清单

1. `pwsh -File tools/check-posts.ps1` 无 ERROR。
2. 文件名日期 == front matter `date` 日期。
3. 图片全部落在 `assets/images/<文章同名目录>/`，链接无括号、无空格。
4. 分类在允许清单内，标签 3-5 个且无全角逗号。
5. 预览：`bundle exec jekyll s`（Ruby 3.3 + `bundle install`），或推到分支后看 Actions。

## 已知待办

- `_drafts/2026-03-01-计算遗传率LDSC.md`：内容完整但作者写的是「吴玥冉」，正文第三、四节（数据准备）疑似重复；发布前需确认作者署名、日期并去重。
- `_drafts/基因组遗传率计算.md`：只有标题和编辑器占位文本，待重写或合并进 LDSC 那篇。
- `assets/images/加权线性回归/`：与 `assets/images/2024-12-05-加权线性回归/` 重复的孤儿目录，无任何引用，确认后可删。
