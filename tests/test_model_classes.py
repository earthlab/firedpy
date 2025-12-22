import os
import shutil

from mock import patch, Mock
from unittest import TestCase

import numpy as np
import xarray as xr

from firedpy.model_classes import EventGrid, EventPerimeter

PROJECT_DIR = os.path.dirname(os.path.dirname(__file__))


class TestEventGridInitMethodTestCase(TestCase):

    def setUp(self):
        # Creating some dummy data to use in tests
        self.out_dir = os.path.join(PROJECT_DIR, "tests", "test_data_dir")
        self.spatial_param = 5
        self.temporal_param = 11
        self.area_unit = "km^2"
        self.time_unit = "days since 1970-01-01"

        self.event1 = EventPerimeter(1, {(1, 2), (2, 3), (3, 4)})
        self.event2 = EventPerimeter(2, {(4, 5), (5, 6), (6, 7)})

    @patch("src.model_classes.xr.open_dataset")
    def test_init_sets_parameters(self, mock_open):
        # Given
        expected_spatial_param = self.spatial_param
        expected_temporal_param = self.temporal_param
        expected_area_unit = self.area_unit
        expected_time_unit = self.time_unit

        mock_array = Mock()
        mock_array.coords = None
        mock_array.value = None

        mock_open.return_value = mock_array

        # When
        grid = EventGrid(self.out_dir, self.spatial_param, self.temporal_param,
                         self.area_unit, self.time_unit)

        # Then
        self.assertEqual(grid._spatial_param, expected_spatial_param)
        self.assertEqual(grid._temporal_param, expected_temporal_param)
        self.assertEqual(grid._area_unit, expected_area_unit)
        self.assertEqual(grid._time_unit, expected_time_unit)

    @patch("src.model_classes.xr.open_dataset")
    def test_init_sets_default_parameters(self, mock_open):
        mock_array = Mock()
        mock_array.coords = None
        mock_array.value = None

        mock_open.return_value = mock_array

        # When
        grid = EventGrid(self.out_dir)

        # Then
        self.assertEqual(grid._spatial_param, 5)
        self.assertEqual(grid._temporal_param, 11)
        self.assertEqual(grid._area_unit, "Unknown")
        self.assertEqual(grid._time_unit, "days since 1970-01-01")

    @patch("src.model_classes.xr.open_dataset")
    def test_merge_perimeters(self, mock_open):
        mock_array = Mock()
        mock_array.coords = None
        mock_array.value = None

        mock_open.return_value = mock_array

        grid = EventGrid(self.out_dir)

        # Store initial coords to compare later
        initial_coords_event1 = self.event1.get_coords().copy()
        initial_coords_event2 = self.event2.get_coords().copy()

        # Method `merge_perimeters`should merge event2 into event1
        merged_perimeters = grid._merge_perimeters(
            [self.event1, self.event2],
            self.event1.get_event_id(),
            self.event2.get_event_id()
        )

        # Validate that the merged perimeter has all the points
        self.assertEqual(
            merged_perimeters[0].get_coords(),
            initial_coords_event1 + initial_coords_event2,
            "Merged coordinates are incorrect"
        )

        # Validate the obsolete event's merge_id and coords
        self.assertEqual(
            merged_perimeters[1].get_merge_id(),
            self.event1.get_event_id(),
            "Merge ID in obsolete event is not set correctly"
        )
        self.assertEqual(
            merged_perimeters[1].get_coords()[0],
            "Merged with event {}".format(self.event1.get_event_id()),
            "Obsolete event coords do not indicate a merge"
        )
        self.assertEqual(
            merged_perimeters[1].get_coords()[1],
            initial_coords_event2,
            "Obsolete event coords do not show original coords"
        )

    @patch("src.model_classes.xr.open_dataset")
    def test_get_available_cells(self, mock_open):
        mock_array = Mock()
        mock_array.coords = None
        mock_array.value = None

        mock_open.return_value = xr.DataArray(
            np.array([[[0, 0], [1, 2]], [[0, 0], [1, 1]]]),
            dims=["time", "y", "x"]
        )

        grid = EventGrid(self.out_dir)

        # Check _get_available_cells method
        available_cells = grid._get_available_cells()

        expected_cells = [(1, 0), (1, 1)]
        self.assertEqual(
            available_cells,
            expected_cells,
            "Available cells do not match expected cells"
        )

    @patch("src.model_classes.xr.open_dataset")
    def test_get_perimeters(self, mock_open):
        mock_open.return_value = xr.DataArray(
            np.array([
                [[0., 0., 0., 0., 0., 0.],
                 [0., 1., 1., 0., 0., 0.],
                 [0., 1., 0., 0., 0., 0.],
                 [0., 0., 0., 0., 0., 0.],
                 [0., 0., 0., 0., 0., 0.],
                 [0., 0., 0., 0., 0., 0.]],

                [[0., 0., 0., 0., 0., 0.],
                 [0., 2., 2., 2., 0., 0.],
                 [0., 2., 2., 0., 0., 0.],
                 [0., 2., 0., 0., 0., 0.],
                 [0., 0., 0., 0., 0., 0.],
                 [0., 0., 0., 0., 0., 0.]],

                [[0., 0., 0., 0., 0., 0.],
                 [0., 0., 3., 3., 0., 0.],
                 [0., 3., 3., 3., 0., 0.],
                 [0., 3., 3., 3., 0., 0.],
                 [0., 3., 3., 3., 0., 0.],
                 [0., 3., 3., 3., 0., 0.]],

                [[0., 0., 0., 0., 0., 0.],
                 [0., 0., 90., 90., 0., 0.],
                 [0., 90., 90., 90., 0., 0.],
                 [0., 90., 90., 90., 0., 0.],
                 [0., 90., 90., 90., 0., 0.],
                 [0., 90., 90., 90., 0., 0.]],

            ]),
            dims=["time", "y", "x"]
        )

        grid = EventGrid(self.out_dir)

        # Check _get_available_cells method
        perimeters = grid.get_event_perimeters()

        self.assertEqual(2, len(perimeters))

        print(perimeters[0].spacetime_coordinates)

    @patch("src.model_classes.xr.open_dataset")
    def test_get_no_available_cells(self, mock_open):
        mock_array = Mock()
        mock_array.coords = None
        mock_array.value = None

        mock_open.return_value = xr.DataArray(
            np.array([[[0, 0], [0, 0]], [[0, 0], [0, 0]]]),
            dims=["time", "y", "x"]
        )

        grid = EventGrid(self.out_dir)

        # Check _get_available_cells method
        available_cells = grid._get_available_cells()

        expected_cells = []
        self.assertEqual(
            available_cells,
            expected_cells,
            "Available cells do not match expected cells"
        )

    @patch("src.model_classes.xr.open_dataset")
    def test_get_all_available_cells(self, mock_open):
        mock_array = Mock()
        mock_array.coords = None
        mock_array.value = None

        mock_open.return_value = xr.DataArray(
            np.array([[[1, 1], [1, 1]], [[1, 1], [1, 1]]]),
            dims=["time", "y", "x"]
        )

        grid = EventGrid(self.out_dir)

        # Check _get_available_cells method
        available_cells = grid._get_available_cells()

        expected_cells = [(0, 0), (0, 1), (1, 0), (1, 1)]
        self.assertEqual(
            available_cells,
            expected_cells,
            "Available cells do not match expected cells"
        )


class TestEventGridSpatialWindow(TestCase):
    _test_data_dir = os.path.join(PROJECT_DIR, "tests", "test_data_dir")

    def tearDown(self) -> None:
        if os.path.exists(self._test_data_dir):
            shutil.rmtree(self._test_data_dir)

    @patch("src.model_classes.xr.open_dataset")
    def test_get_spatial_window_center(self, mock_open):
        # Assuming we have an array of 11x11 (y, x)
        array_dims = (11, 11)

        mock_array = Mock()
        mock_array.coords = None
        mock_array.value = None

        mock_open.return_value = mock_array

        # Case where the window is centered (i.e., not on edge of the grid.)
        y, x = (5, 5)
        expected_output = (3, 7, 3, 7, [2, 2], [3, 3])  # Based on event logic

        event_grid = EventGrid(self._test_data_dir, 2)

        result = event_grid._get_spatial_window(y, x, array_dims)
        msg = (f"For (y, x) = ({y}, {x}), expected {expected_output} but got "
               f"{result}")
        self.assertEqual(result, expected_output, msg)

    @patch("src.model_classes.xr.open_dataset")
    def test_get_spatial_window_top_left(self, mock_open):
        # Testing a case where the window is in the top-left corner.
        array_dims = (11, 11)

        mock_array = Mock()
        mock_array.coords = None
        mock_array.value = None

        mock_open.return_value = mock_array

        y, x = (1, 1)
        expected_output = (0, 3, 0, 3, [1, 1], [0, 0])

        event_grid = EventGrid(self._test_data_dir, 2)

        result = event_grid._get_spatial_window(y, x, array_dims)
        msg = (f"For (y, x) = ({y}, {x}), expected {expected_output} but got "
               f"{result}")
        self.assertEqual(result, expected_output, msg)

    @patch("src.model_classes.xr.open_dataset")
    def test_get_spatial_window_top_right(self, mock_open):
        # Testing a case where the window is in the top-left corner.
        array_dims = (11, 11)

        mock_array = Mock()
        mock_array.coords = None
        mock_array.value = None

        mock_open.return_value = mock_array

        y, x = (0, 10)
        expected_output = (0, 2, 8, 11, [0, 2], [0, 8])

        event_grid = EventGrid(self._test_data_dir, 2)

        result = event_grid._get_spatial_window(y, x, array_dims)
        msg = (f"For (y, x) = ({y}, {x}), expected {expected_output} but got "
               f"{result}")
        self.assertEqual(result, expected_output, msg)

    @patch("src.model_classes.xr.open_dataset")
    def test_get_spatial_window_bottom_right(self, mock_open):
        # Testing a case where the window is in the bottom-right corner.
        array_dims = (11, 11)

        mock_array = Mock()
        mock_array.coords = None
        mock_array.value = None

        mock_open.return_value = mock_array

        y, x = (10, 10)
        expected_output = (8, 11, 8, 11, [2, 2], [8, 8])

        event_grid = EventGrid(self._test_data_dir, 2)

        result = event_grid._get_spatial_window(y, x, array_dims)
        msg = (f"For (y, x) = ({y}, {x}), expected {expected_output} but got "
               f"{result}")
        self.assertEqual(result, expected_output, msg)

    @patch("src.model_classes.xr.open_dataset")
    def test_get_spatial_window_bottom_left(self, mock_open):
        # Testing a case where the window is in the top-left corner.
        array_dims = (11, 11)

        mock_array = Mock()
        mock_array.coords = None
        mock_array.value = None

        mock_open.return_value = mock_array

        y, x = (10, 1)
        expected_output = (8, 11, 0, 3, [2, 1], [8, 0])

        event_grid = EventGrid(self._test_data_dir, 2)

        result = event_grid._get_spatial_window(y, x, array_dims)
        msg = (f"For (y, x) = ({y}, {x}), expected {expected_output} but got "
               f"{result}")
        self.assertEqual(result, expected_output, msg)

    @patch("src.model_classes.xr.open_dataset")
    def test_get_spatial_window_top_left_corner(self, mock_open):
        # Testing a case where the window is in the top-left corner.
        array_dims = (11, 11)

        mock_array = Mock()
        mock_array.coords = None
        mock_array.value = None

        mock_open.return_value = mock_array

        y, x = (0, 0)
        expected_output = (0, 2, 0, 2, [0, 0], [0, 0])

        event_grid = EventGrid(self._test_data_dir, 2)

        result = event_grid._get_spatial_window(y, x, array_dims)
        msg = (f"For (y, x) = ({y}, {x}), expected {expected_output} but got "
               f"{result}")
        self.assertEqual(result, expected_output, msg)

    @patch("src.model_classes.xr.open_dataset")
    def test_get_spatial_window_top_right_corner(self, mock_open):
        # Testing a case where the window is in the top-left corner.
        array_dims = (11, 11)

        mock_array = Mock()
        mock_array.coords = None
        mock_array.value = None

        mock_open.return_value = mock_array

        y, x = (0, 11)
        expected_output = (0, 2, 9, 11, [0, 2], [0, 9])

        event_grid = EventGrid(self._test_data_dir, 2)

        result = event_grid._get_spatial_window(y, x, array_dims)
        msg = (f"For (y, x) = ({y}, {x}), expected {expected_output} but got "
               f"{result}")
        self.assertEqual(result, expected_output, msg)

    @patch("src.model_classes.xr.open_dataset")
    def test_get_spatial_window_bottom_right_corner(self, mock_open):
        # Testing a case where the window is in the bottom-right corner.
        array_dims = (11, 11)

        mock_array = Mock()
        mock_array.coords = None
        mock_array.value = None

        mock_open.return_value = mock_array

        y, x = (11, 11)
        expected_output = (9, 11, 9, 11, [2, 2], [9, 9])

        event_grid = EventGrid(self._test_data_dir, 2)

        result = event_grid._get_spatial_window(y, x, array_dims)
        msg = (f"For (y, x) = ({y}, {x}), expected {expected_output} but got "
               f"{result}")
        self.assertEqual(result, expected_output, msg)

    @patch("src.model_classes.xr.open_dataset")
    def test_get_spatial_window_bottom_left_corner(self, mock_open):
        # Testing a case where the window is in the top-left corner.
        array_dims = (11, 11)

        mock_array = Mock()
        mock_array.coords = None
        mock_array.value = None

        mock_open.return_value = mock_array

        y, x = (11, 0)
        expected_output = (9, 11, 0, 2, [2, 0], [9, 0])

        event_grid = EventGrid(self._test_data_dir, 2)

        result = event_grid._get_spatial_window(y, x, array_dims)
        msg = (f"For (y, x) = ({y}, {x}), expected {expected_output} but got "
               f"{result}")
        self.assertEqual(result, expected_output, msg)
