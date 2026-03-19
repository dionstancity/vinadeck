from __future__ import annotations

import json
from uuid import uuid4

import streamlit.components.v1 as components

from dockgui.models import DockingBox


def render_structure_viewer(
    receptor_pdb: str | None,
    ligand_pdb: str | None,
    docking_box: DockingBox | None,
    *,
    height: int = 560,
) -> None:
    viewer_id = f"viewer-{uuid4().hex}"
    receptor_json = json.dumps(receptor_pdb or "")
    ligand_json = json.dumps(ligand_pdb or "")
    box_json = json.dumps(docking_box.as_dict() if docking_box else None)

    html = f"""
    <div id="{viewer_id}" style="width: 100%; height: {height}px; position: relative;"></div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
      const receptor = {receptor_json};
      const ligand = {ligand_json};
      const box = {box_json};
      const viewer = $3Dmol.createViewer(document.getElementById("{viewer_id}"), {{ backgroundColor: "white" }});

      if (receptor) {{
        viewer.addModel(receptor, "pdb");
        viewer.setStyle({{model: 0}}, {{cartoon: {{color: "spectrum"}}, line: {{hidden: true}}}});
      }}

      if (ligand) {{
        viewer.addModel(ligand, "pdb");
        const ligandModelIndex = receptor ? 1 : 0;
        viewer.setStyle({{model: ligandModelIndex}}, {{stick: {{radius: 0.18, colorscheme: "cyanCarbon"}}, sphere: {{radius: 0.28, colorscheme: "cyanCarbon"}}}});
      }}

      if (box) {{
        viewer.addBox({{
          center: {{x: box.center_x, y: box.center_y, z: box.center_z}},
          dimensions: {{w: box.size_x, h: box.size_y, d: box.size_z}},
          color: "magenta",
          alpha: 0.2,
          wireframe: true
        }});
      }}

      viewer.zoomTo();
      viewer.render();
    </script>
    """
    components.html(html, height=height)

