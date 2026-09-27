# -*- coding: utf-8 -*-

"""Tests for where PEMT stores its data."""

import os

import pemt.constants as constants


def test_env_variable_wins(monkeypatch, tmp_path):
    monkeypatch.setenv("PEMT_DATA_DIR", str(tmp_path / "custom"))
    assert constants.default_data_dir() == str(tmp_path / "custom")


def test_source_checkout_uses_repo_data(monkeypatch):
    monkeypatch.delenv("PEMT_DATA_DIR", raising=False)
    repo_root = os.path.normpath(os.path.join(constants.HERE, "..", ".."))
    assert constants.default_data_dir() == os.path.join(repo_root, "data")


def test_installed_package_uses_working_directory(monkeypatch, tmp_path):
    monkeypatch.delenv("PEMT_DATA_DIR", raising=False)
    fake_install = tmp_path / "site-packages" / "pemt"
    fake_install.mkdir(parents=True)
    monkeypatch.setattr(constants, "HERE", str(fake_install))
    monkeypatch.chdir(tmp_path)
    assert constants.default_data_dir() == str(tmp_path / "pemt_data")
