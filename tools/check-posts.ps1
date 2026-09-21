# 校验 _posts / _drafts 的文章规范：文件名、front matter、日期、分类标签、图片链接
# 用法: pwsh -File tools/check-posts.ps1

$ErrorActionPreference = 'Stop'
$repo = Split-Path -Parent $PSScriptRoot
Set-Location $repo

$script:errCount = 0
$script:warnCount = 0
$script:infoCount = 0

function Write-Fail([string]$Message) {
  $script:errCount++
  Write-Host "ERROR   $Message" -ForegroundColor Red
}

function Write-Warn([string]$Message) {
  $script:warnCount++
  Write-Host "WARN    $Message" -ForegroundColor Yellow
}

function Write-Info([string]$Message) {
  $script:infoCount++
  Write-Host "INFO    $Message" -ForegroundColor DarkGray
}

# 历史遗留：文件名带空格，改名的收益小于旧链接失效的风险，保留但记录在此
$legacyNames = @('2025-01-03-call variants.md', '2025-06-24-docker 简单教程.md')

$allowed = [ordered]@{
  '生物信息学' = @('测序与变异', '基因组分析')
  '统计方法'   = @('回归模型', '蒙特卡洛方法')
  '统计遗传学' = @('孟德尔随机化', '遗传率估计')
  '软件工具'   = @('环境管理')
  '学术写作'   = @('论文各节')
}

function Get-FrontMatter([string]$Path) {
  $lines = [System.IO.File]::ReadAllLines($Path)
  if ($lines.Count -eq 0 -or $lines[0].Trim() -ne '---') { return $null }
  $end = -1
  for ($i = 1; $i -lt $lines.Count; $i++) {
    if ($lines[$i].Trim() -eq '---') { $end = $i; break }
  }
  if ($end -lt 0) { return $null }
  $map = @{}
  for ($i = 1; $i -lt $end; $i++) {
    $line = $lines[$i]
    if ($line.Trim() -eq '' -or $line -match '^\s*#') { continue }
    if ($line -match '^([A-Za-z_][A-Za-z0-9_]*)\s*:\s*(.*)$') {
      $map[$Matches[1]] = $Matches[2].Trim()
    }
  }
  return $map
}

function Parse-List([string]$Value) {
  if ([string]::IsNullOrWhiteSpace($Value)) { return @() }
  $v = $Value.Trim()
  if ($v.StartsWith('[') -and $v.EndsWith(']')) { $v = $v.Substring(1, $v.Length - 2) }
  return @($v -split ',' | ForEach-Object { $_.Trim().Trim('"').Trim("'") } | Where-Object { $_ -ne '' })
}

$postsDir = Join-Path $repo '_posts'
$draftsDir = Join-Path $repo '_drafts'
$posts = @(Get-ChildItem -Path $postsDir -Filter '*.md' -File -ErrorAction SilentlyContinue)
$drafts = @(Get-ChildItem -Path $draftsDir -Filter '*.md' -File -ErrorAction SilentlyContinue)

$tagCount = @{}
$catCount = @{}
$pinned = 0

foreach ($file in $posts) {
  $name = $file.Name
  $rel = "_posts/$name"

  if ($name -notmatch '^(\d{4}-\d{2}-\d{2})-(.+)\.md$') {
    Write-Fail "$rel 缺少 YYYY-MM-DD- 前缀（Jekyll 不会构建这篇，线上看不到）"
    continue
  }
  $fileDate = $Matches[1]
  $slug = $Matches[2]

  $isLegacy = $legacyNames -contains $name

  if ($slug -match '[ ()\u3000：:／/]') {
    if ($isLegacy) {
      Write-Info "$rel 文件名含空格（历史遗留，白名单内）"
    }
    else {
      Write-Warn "$rel 文件名含空格或符号，会进入 permalink /posts/:title/"
    }
  }

  $fm = Get-FrontMatter $file.FullName
  if ($null -eq $fm) {
    Write-Fail "$rel 缺少 front matter（文件需以 --- 开头、以 --- 结束）"
    continue
  }

  foreach ($key in @('title', 'description', 'date', 'categories', 'tags')) {
    if (-not $fm.ContainsKey($key) -or $fm[$key].Trim() -eq '') {
      Write-Fail "$rel front matter 缺少 $key"
    }
  }

  if ($fm.ContainsKey('date') -and $fm['date'] -match '^(\d{4}-\d{2}-\d{2})') {
    if ($Matches[1] -ne $fileDate) {
      Write-Fail "$rel 文件名日期 $fileDate 与 front matter date $($Matches[1]) 不一致"
    }
  }

  $categories = @(Parse-List $fm['categories'])
  if ($categories.Count -ne 2) {
    Write-Fail "$rel categories 应为一级+二级共 2 项，当前 $($categories.Count) 项: [$($categories -join ', ')]"
  }
  elseif (-not $allowed.Contains($categories[0])) {
    Write-Fail "$rel 一级分类 '$($categories[0])' 不在允许清单: $($allowed.Keys -join ' / ')"
  }
  elseif ($allowed[$categories[0]] -notcontains $categories[1]) {
    Write-Fail "$rel 二级分类 '$($categories[1])' 不属于 '$($categories[0])'，允许: $($allowed[$categories[0]] -join ' / ')"
  }
  else {
    $k = "$($categories[0])/$($categories[1])"
    $catCount[$k] = 1 + [int]$catCount[$k]
  }

  $rawTags = [string]$fm['tags']
  if ($rawTags -match '，') {
    Write-Fail "$rel tags 含全角逗号，会被解析成一个假标签"
  }
  $tags = @(Parse-List $rawTags)
  if ($tags.Count -lt 3 -or $tags.Count -gt 5) {
    Write-Warn "$rel 标签 $($tags.Count) 个（建议 3-5）: [$($tags -join ', ')]"
  }
  foreach ($t in $tags) {
    if ($categories -contains $t) { Write-Warn "$rel 标签 '$t' 与分类名重复" }
    $tagCount[$t] = 1 + [int]$tagCount[$t]
  }

  if ($fm.ContainsKey('pin') -and $fm['pin'] -match '^(true|yes)$') {
    $pinned++
  }

  $text = [System.IO.File]::ReadAllText($file.FullName)
  foreach ($m in [regex]::Matches($text, '!\[[^\]]*\]\((.*?)\)\s*$', 'Multiline')) {
    $ref = $m.Groups[1].Value.Trim()
    if ($ref -match '^(https?:|//|/|#|data:)') { continue }
    $clean = [uri]::UnescapeDataString(($ref -replace '^\.\./', ''))
    if (-not (Test-Path -LiteralPath (Join-Path $repo $clean))) {
      Write-Fail "$rel 图片不存在: $ref"
    }
    if (($clean -match '[ ()]') -and (-not $isLegacy)) {
      Write-Warn "$rel 图片路径含空格或括号，建议改名: $ref"
    }
  }
}

foreach ($file in $drafts) {
  $fm = Get-FrontMatter $file.FullName
  if ($null -eq $fm) {
    Write-Warn "_drafts/$($file.Name) 缺少 front matter，发布前需补"
    continue
  }
  if (-not $fm.ContainsKey('date')) {
    Write-Warn "_drafts/$($file.Name) 未写 date，移到 _posts 前需补齐并加文件名日期前缀"
  }
}

Write-Host ''
Write-Host "文章: $($posts.Count) 已发布, $($drafts.Count) 草稿" -ForegroundColor Cyan
Write-Host '分类分布:'
foreach ($k in ($catCount.Keys | Sort-Object)) { Write-Host ("  {0,-26} {1}" -f $k, $catCount[$k]) }
Write-Host '标签使用次数:'
foreach ($k in ($tagCount.Keys | Sort-Object { -$tagCount[$_] })) { Write-Host ("  {0,-26} {1}" -f $k, $tagCount[$k]) }

Write-Host ''
if ($script:errCount -gt 0) {
  Write-Host "FAIL: $($script:errCount) 个错误, $($script:warnCount) 个警告" -ForegroundColor Red
  exit 1
}
Write-Host "PASS: 0 个错误, $($script:warnCount) 个警告" -ForegroundColor Green
