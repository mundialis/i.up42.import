#!/usr/bin/env python3
#
############################################################################
#
# MODULE:       i.up42.import
# AUTHOR(S):    Anika Weinmann
# PURPOSE:      Imports Pléiades data using the Python SDK from UP42
#
# COPYRIGHT:    (C) 2020-2022 by mundialis GmbH & Co. KG and the GRASS Development Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
############################################################################

# %module
# % description: Imports Pléiades data using the Python SDK from UP42.
# % keyword: imagery
# % keyword: import
# % keyword: projection
# % keyword: Pléiades
# % keyword: UP42
# % keyword: VHR
# %end

# %option G_OPT_V_INPUT
# % key: input
# % required: no
# % description: Input vector map of BBOX
# %end

# %option
# % key: input_geojson
# % type: string
# % required: no
# % multiple: no
# % description: Input map as GeoJSON string in EPSGS:4326
# %end

# %option
# % key: producttype
# % type: string
# % required: no
# % multiple: no
# % options: pleiades
# % description: UP42 product type to filter
# % answer: pleiades
# % guisection: Filter
# %end

# %option
# % key: start
# % type: string
# % required: no
# % multiple: no
# % description: Start date ('YYYY-MM-DD')
# % guisection: Filter
# % answer: 2020-01-01
# %end

# %option
# % key: end
# % type: string
# % required: no
# % multiple: no
# % description: End date ('YYYY-MM-DD')
# % guisection: Filter
# % answer: today
# %end

# %option
# % key: clouds
# % type: integer
# % required: no
# % multiple: no
# % description: Maximum cloud cover percentage for scene
# % guisection: Filter
# % answer: 100
# %end

# %option
# % key: area_relation
# % type: string
# % required: no
# % multiple: no
# % options: Intersects,Contains
# % description: Spatial relation
# % answer: Contains
# % guisection: Region
# %end

# %option
# % key: directory
# % type: string
# % required: no
# % multiple: no
# % description: Directory to save the original data; If not set, downloaded data will be deleted
# %end

# %option
# % key: output
# % type: string
# % required: no
# % multiple: no
# % key_desc: name
# % description: Name for output raster map
# % gisprompt: new,cell,raster
# % guisection: Output
# %end

# %option
# % key: max_area
# % type: integer
# % required: no
# % multiple: no
# % description: Maximum allowed area size in sqm
# % guisection: Filter
# %end

# %flag
# % key: p
# % description: Print only credits for job
# %end

# %flag
# % key: m
# % description: Create vector map with meta data named meta_{output}
# %end

# %rules
# % required: input,input_geojson
# % requires: -m, output
# %end


import atexit
from datetime import date
import json
import os
import sys
import up42

import grass.script as grass

# initialize global vars
tmpfolder = None
rm_vectors = []
rm_rasters = []


def cleanup():
    grass.message(_("Cleaning up..."))
    nuldev = open(os.devnull, "w")
    kwargs = {"flags": "f", "quiet": True, "stderr": nuldev}
    for rmrast in rm_rasters:
        if grass.find_file(name=rmrast, element="raster")["file"]:
            grass.run_command("g.remove", type="raster", name=rmrast, **kwargs)
    for rmvect in rm_vectors:
        if grass.find_file(name=rmvect, element="vector")["file"]:
            grass.run_command("g.remove", type="vector", name=rmvect, **kwargs)
    if tmpfolder:
        grass.try_rmdir(os.path.join(tmpfolder))


def check_start_end(start, end):
    # set end if end is today
    if end == "today":
        end = date.today().strftime("%Y-%m-%d")
    # check if end is after start
    start_list = [int(x) for x in start.split("-")]
    start_date = date(start_list[0], start_list[1], start_list[2])
    end_list = [int(x) for x in end.split("-")]
    end_date = date(end_list[0], end_list[1], end_list[2])
    if end_date < start_date:
        grass.fatal(_("End date is before start date"))

    return start_date, end_date


def main():

    global tmpfolder, rm_vectors, rm_rasters

    if options["input"] == "None":
        options["input"] = None
    if options["input_geojson"] == "None":
        options["input_geojson"] = None
    elif "'" in options["input_geojson"]:
        options["input_geojson"] = options["input_geojson"].replace("'", '"')

    # input = options['input']
    producttype = options["producttype"]
    start = options["start"]
    end = options["end"]
    clouds = int(options["clouds"])
    area_relation = options["area_relation"]
    output = options["output"]
    if options["directory"]:
        folder = options["directory"]
        if not os.path.isdir(folder):
            os.makedirs(folder)
    else:
        folder = grass.tempdir()
        tmpfolder = folder

    start_date, end_date = check_start_end(start, end)

    # set some common environmental variables, like:
    os.environ.update(
        dict(
            GRASS_COMPRESS_NULLS="1",
            GRASS_COMPRESSOR="LZ4",
            GRASS_MESSAGE_FORMAT="plain",
        )
    )

    if flags["m"] and options["input_geojson"]:
        if not grass.find_program("v.in.geojson", "--help"):
            grass.message(
                _(
                    "The 'v.in.geojson' module was not found, you can install it with:"
                    + "\ng.extension v.in.geojson url=...\n\n"
                )
            )
    elif options["input"]:
        if not grass.find_program("v.out.geojson", "--help"):
            grass.message(
                _(
                    "The 'v.out.geojson' module was not found, you can install it with:"
                    + "\ng.extension v.out.geojson url=...\n\n"
                )
            )

    # check if max_area is exceeded
    if options["max_area"]:
        maximum = float(options["max_area"])
        tmp_aoi = f"tmp_aoi_{str(os.getpid())}"
        rm_vectors.append(tmp_aoi)
        if options["input_geojson"]:
            grass.run_command(
                "v.in.geojson",
                input=options["input_geojson"],
                output=tmp_aoi,
                quiet=True,
            )
        else:
            grass.run_command("g.copy", vector=f"{options['input']},{tmp_aoi}")
        grass.run_command(
            "v.db.addcolumn",
            map=tmp_aoi,
            columns="tmparea double precision",
            quiet=True,
            overwrite=True,
        )
        grass.run_command(
            "v.to.db",
            map=tmp_aoi,
            option="area",
            columns="tmparea",
            units="meters",
            quiet=True,
            overwrite=True,
        )
        area_sqm = sum(
            [
                float(x)
                for x in grass.parse_command(
                    "v.db.select", map=tmp_aoi, columns="tmparea", flags="c"
                )
            ]
        )
        if area_sqm > maximum:
            grass.fatal(
                _(
                    "The input vector has with %.2f sqm a larger area "
                    "than the given maximum (%.2f sqm)" % (area_sqm, maximum)
                )
            )
        else:
            grass.message(_("The input vector has an area of %.2f sqm") % area_sqm)

    grass.message(_("Authenticate & access UP42 project ..."))
    up42.authenticate(
        project_id=os.environ["UP42_PROJECT_ID"],
        project_api_key=os.environ["UP42_PROJECT_API_KEY"],
    )
    project = up42.initialize_project()

    workflow = project.create_workflow(name="dl_tiff_cloud", use_existing=True)
    # add the workflow tasks
    w_tasks = [
        "oneatlas-pleiades-fullscene",
        "data-conversion-dimap",
        "oneatlas-cloudmask",
    ]
    workflow.add_workflow_tasks(w_tasks)

    # input vector map to GeoJSON and read geometry
    if options["input"]:
        gj = json.loads(
            [
                key
                for key in grass.parse_command(
                    "v.out.geojson", input=options["input"], output="-"
                )
            ][0]
        )
        geom = gj["features"][0]["geometry"]
    elif options["input_geojson"]:
        gj = json.loads(options["input_geojson"])
        geom = gj["features"][0]["geometry"]
        # check if GeoJSON coordinates are in lat lon
        lat = [x[0] for x in geom["coordinates"][0]]
        lon = [x[1] for x in geom["coordinates"][0]]
        if not all((x >= -180 and x <= 180) for x in lat):
            grass.fatal(
                _("Latitude values of GeoJSON are not in range of -180 and 180")
            )
        if not all((x >= -90 and x <= 90) for x in lon):
            grass.fatal(_("Longitude values of GeoJSON are not in range of -90 and 90"))
    else:
        grass.fatal(_("<input> or <input_geojson have to be set>"))

    # set time
    time_str = "%sT00:00:00+00:00/%sT23:59:59+00:00" % (
        start_date.strftime("%Y-%m-%d"),
        end_date.strftime("%Y-%m-%d"),
    )

    if producttype == "pleiades":
        grass.message(_("Pleiades download parameter ..."))
        # https://docs.up42.com/up42-blocks/data/pleiades-reflectance-download.html
        input_parameters = {
            "oneatlas-pleiades-fullscene:1": {
                "ids": None,
                "time": time_str,
                "limit": 1,
                "asset_ids": None,
                "time_series": None,
                "max_cloud_cover": clouds,
                area_relation.lower(): geom,
            },
            "data-conversion-dimap:1": {
                "ms": True,
                "pan": False,
                "bbox": None,
                "contains": None,
                "intersects": None,
                "clip_to_aoi": False,
            },
            "oneatlas-cloudmask:1": {},
        }
    else:
        grass.fatal(_("The producttype <%s> is not supported yet") % producttype)

    # Test parameters
    grass.message(_("Test parameter..."))
    try:
        test_job = workflow.test_job(
            input_parameters=input_parameters, track_status=True
        )
    except Exception as e:
        grass.fatal(_("Error in testing UP42 job: %s") % e)
    if test_job.info["status"] != "SUCCEEDED":
        grass.fatal(
            _("Error in testing UP42 job. Status is <%s>") % test_job.info["status"]
        )

    test_results = test_job.get_results_json()
    for feat in test_results["features"]:
        grass.message(
            _("Found scene with ID <%s> - cloud coverage: %f - from <%s>")
            % (
                feat["properties"]["id"],
                feat["properties"]["cloudCover"],
                feat["properties"]["acquisitionDate"],
            )
        )
    # This prints the estimated credits automatically
    workflow.estimate_job(input_parameters=input_parameters)
    if flags["m"] and options["output"]:
        grass.message(_("Write metadata to vector map..."))
        voutput = "meta_%s" % output
        if options["input_geojson"]:
            grass.run_command(
                "v.in.geojson",
                input=options["input_geojson"],
                output=voutput,
                quiet=True,
            )
        else:
            grass.run_command("g.copy", vector="%s,%s" % (options["input"], voutput))
        if len(test_results["features"]) > 1:
            grass.fatal(
                _("Metadata vector map only implemented for one resulting feature")
            )
        metadata = test_results["features"][0]["properties"]
        addcolumn = []
        for key, val in metadata.items():
            if isinstance(val, str):
                addcolumn.append("%s VARCHAR(15)" % key)
            elif isinstance(val, float):
                addcolumn.append("%s DOUBLE PRECISION" % key)
            elif isinstance(val, int):
                addcolumn.append("%s INT" % key)
        if len(grass.vector_db(voutput)) == 0:
            grass.run_command("v.db.addtable", map=voutput, quiet=True)
        grass.run_command(
            "v.db.addcolumn", map=voutput, columns=",".join(addcolumn), quiet=True
        )
        for key, val in metadata.items():
            newval = val
            # bool area also int but v,db.update gives error
            if isinstance(val, int):
                newval = int(val)
            grass.run_command(
                "v.db.update",
                map=voutput,
                layer=1,
                column=key,
                value=newval,
                quiet=True,
            )
        grass.message(_("Metadata to vector map <%s> created" % voutput))

    # Download scene
    if not flags["p"]:
        grass.message(_("Running job..."))
        job = workflow.run_job(input_parameters=input_parameters, track_status=True)
        grass.message(_("Downloading job result..."))
        jobid = job.info["id"]
        jobtasks = job.get_jobtasks()[0].info
        jobtasksids = [tsk["id"] for tsk in jobtasks]

        dimapfolder = os.path.join(folder, "dimap")
        tiffolder = os.path.join(folder, "tif")
        for fol in [dimapfolder, tiffolder]:
            try:
                if not os.path.isdir(fol):
                    os.makedirs(fol)
            except Exception as e:
                grass.fatal(_("Could not create directory {}:" " {}".format(fol, e)))
        # in the DIMAP--> TIFF conversion, the PAN is lost, so the first
        # task result has to be downloaded individually
        # download the MS + PAN in DIMAP format (result of first task):
        jobtask1 = up42.initialize_jobtask(jobtask_id=jobtasksids[0], job_id=jobid)
        results = jobtask1.download_results(dimapfolder)
        for result in results:
            name = os.path.basename(result)
            if name.startswith("DIM") and name.endswith(".XML"):
                if "_MS_" in name:
                    grass.run_command("r.import", input=result, output=output)
                    grass.run_command(
                        "g.rename", raster="{}.1,{}.red".format(output, output)
                    )
                    grass.run_command(
                        "g.rename", raster="{}.2,{}.green".format(output, output)
                    )
                    grass.run_command(
                        "g.rename", raster="{}.3,{}.blue".format(output, output)
                    )
                    grass.run_command(
                        "g.rename", raster="{}.4,{}.nir".format(output, output)
                    )
                elif "_P_" in name:
                    grass.run_command(
                        "r.import", input=result, output="{}.pan".format(output)
                    )

        # download final results in tif; here only the cloud mask is needed
        results_tif = job.download_results(tiffolder)
        rastername = "tmp_import_raster_{}".format(os.getpid())
        for result in results_tif:
            if result.endswith(".tif"):
                grass.run_command("r.import", input=result, output=rastername)
                grass.run_command(
                    "g.rename", raster="{}.5,{}.alpha".format(rastername, output)
                )
                for i in range(1, 5):
                    rm_rasters.append("{}.{}".format(rastername, i))

        grass.run_command(
            "i.group",
            group=output,
            input=("%s.red,%s.green," "%s.blue,%s.nir,%s.pan,%s.alpha")
            % (output, output, output, output, output, output),
        )
        grass.message(_("Created group <%s>") % output)

    return 0


if __name__ == "__main__":
    options, flags = grass.parser()
    atexit.register(cleanup)
    sys.exit(main())
