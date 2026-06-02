"""Tool functions for the pdf sub-skill."""
from __future__ import annotations

from .plans.T1_font_embed import pdf_font_embed_check  # noqa: F401
from .plans.T2_inspect import pdf_inspect  # noqa: F401
from .plans.T3_merge_split import pdf_merge, pdf_split, pdf_extract_pages  # noqa: F401
from .plans.T4_crop_whitespace import pdf_crop_whitespace  # noqa: F401
from .plans.T5_bookmarks import pdf_list_bookmarks, pdf_set_bookmarks  # noqa: F401
from .plans.T6_metadata import pdf_get_metadata, pdf_set_metadata  # noqa: F401
from .plans.T7_watermark import pdf_watermark_add  # noqa: F401
from .plans.T8_text_extract import pdf_extract_text  # noqa: F401
from .plans.T9_compress import pdf_compress  # noqa: F401
from .plans.T10_health_report import pdf_health_report  # noqa: F401
from .plans.T13_lint_dpi import pdf_lint_image_detection_dpi  # noqa: F401
