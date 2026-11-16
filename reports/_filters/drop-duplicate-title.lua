--[[
  Drop a body H1 that merely repeats the page title.

  Every report page here is a thin wrapper: the .qmd sets `title:` in its front
  matter and then `{{< include >}}`s an analysis REPORT.md that opens with the
  same title as a level-1 heading, because that heading is what makes the file
  readable on GitHub. Rendered, the title therefore appeared twice -- once in
  Quarto's title block and again at the top of the body.

  Fixing it in REPORT.md would break the GitHub view, and fixing it by dropping
  `title:` would cost the title block, the navbar text and the <title> element.
  So it is fixed here instead.

  The guard is deliberately narrow: only the FIRST level-1 heading, and only
  when its text is exactly the document title. A report whose first H1 says
  something else keeps it, and H1s further down are never touched.
]]

function Pandoc(doc)
  if not doc.meta.title then return doc end
  local title = pandoc.utils.stringify(doc.meta.title)

  for i, blk in ipairs(doc.blocks) do
    if blk.t == "Header" and blk.level == 1 then
      if pandoc.utils.stringify(blk.content) == title then
        table.remove(doc.blocks, i)
      end
      break -- only ever consider the first H1
    end
  end
  return doc
end
